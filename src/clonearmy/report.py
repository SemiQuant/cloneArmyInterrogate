from typing import Dict, List, Set, Tuple
import json
from pathlib import Path
import base64
from io import BytesIO, StringIO
import datetime
import logging
from collections import defaultdict
import shutil
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import click

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from jinja2 import Template
import os

# Configure plotting styles properly
sns.set_theme(style="whitegrid")  # This sets up seaborn's styling
plt.style.use('default')  # Use matplotlib's default style as a base

# Disable interactive mode
plt.ioff()

from .processor import AmpliconProcessor

logger = logging.getLogger(__name__)

def get_mutations_from_haplotype(haplotype: str, reference_seq: str, count: int,
                                 insertions: str = '') -> List[dict]:
    """Extract mutations (including insertions) from a haplotype sequence."""
    mutations = []
    for pos, (ref, var) in enumerate(zip(reference_seq.upper(), haplotype)):
        if var.islower() or var == '-':  # Mutation or deletion
            mutations.append({
                'position': pos + 1,
                'reference_base': ref,
                'mutation': var.upper(),
                'count': count
            })
    for pos, bases in AmpliconProcessor.parse_insertions(insertions):
        mutations.append({
            'position': pos,
            'reference_base': '-',
            'mutation': f"+{bases}",
            'count': count
        })
    return mutations

def fig_to_base64() -> str:
    """Convert matplotlib figure to base64 string."""
    try:
        buf = BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', dpi=300, backend='Agg')
        buf.seek(0)
        plt.close()
        return f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}"
    except Exception as e:
        logger.error(f"Error converting figure to base64: {str(e)}")
        return ""

def create_mutation_frequency_plot(results_df: pd.DataFrame, reference_seq: str) -> str:
    """Create an interactive mutation frequency plot using Plotly."""
    if results_df.empty:
        return ""
        
    # Process mutations from all haplotypes
    all_mutations = []
    total_reads = results_df['count'].sum()
    
    if total_reads == 0:
        return ""
    
    for _, row in results_df.iterrows():
        mutations = get_mutations_from_haplotype(row['haplotype'], reference_seq, row['count'],
                                                 row.get('insertions', ''))
        all_mutations.extend(mutations)
    
    if not all_mutations:
        return ""
    
    # Aggregate mutation frequencies
    mutation_freq = defaultdict(int)
    for mut in all_mutations:
        key = f"{mut['position']} {mut['reference_base']}→{mut['mutation']}"
        mutation_freq[key] += mut['count']
    
    if not mutation_freq:
        return ""
    
    # Convert to percentages and sort by position
    positions = []
    frequencies = []
    labels = []
    
    for key, count in sorted(mutation_freq.items(), key=lambda x: int(x[0].split()[0])):
        positions.append(key)
        freq = (count / total_reads) * 100
        frequencies.append(freq)
        labels.append(f"{key}: {freq:.2f}%")
    
    if not frequencies:
        return ""
    
    # Create Plotly figure
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=positions,
        y=frequencies,
        marker_color='rgb(31, 119, 180)',
        hovertext=labels,
        name='Mutation Frequency'
    ))
    
    max_freq = max(frequencies) if frequencies else 0
    
    fig.update_layout(
        title='Mutation Frequencies',
        xaxis=dict(
            title='Position and Mutation',
            tickangle=45,
            tickfont=dict(size=10)
        ),
        yaxis=dict(
            title='Frequency (%)',
            range=[0, max_freq * 1.1] if max_freq > 0 else [0, 100]
        ),
        margin=dict(b=100, l=60, r=20, t=40),
        width=max(800, len(positions) * 50) if positions else 800,
        height=500,
        template='plotly_white'
    )
    
    # Convert to HTML
    plot_html = fig.to_html(
        full_html=False,
        include_plotlyjs='cdn',
        config={'displayModeBar': True}
    )
    
    return plot_html

def create_mutation_spectrum(results: Dict[str, pd.DataFrame], reference_seq: str) -> str:
    """Create mutation spectrum analysis."""
    try:
        mutation_types = defaultdict(int)
        total_mutations = 0
        double_mutations = defaultdict(int)
        
        # Normalize reference sequence once
        reference_seq = reference_seq.upper()
        
        for df in results.values():
            if df.empty:
                continue
                
            for _, row in df.iterrows():
                haplotype = row['haplotype']
                count = row['count']
                
                # Track positions with mutations for linked analysis
                mutation_positions = []
                
                for i, (ref, var) in enumerate(zip(reference_seq, haplotype)):
                    if var.islower():  # This identifies a mutation
                        # Create mutation string with both bases in uppercase
                        mutation = f"{ref}>{var.upper()}"
                        # Only count if they're actually different after normalization
                        if ref != var.upper():
                            mutation_types[mutation] += count
                            total_mutations += count
                            mutation_positions.append((i + 1, ref, var.upper()))
                
                # Process linked mutations
                for i in range(len(mutation_positions)):
                    for j in range(i + 1, len(mutation_positions)):
                        pos1, ref1, mut1 = mutation_positions[i]
                        pos2, ref2, mut2 = mutation_positions[j]
                        double_key = f"({pos1}{ref1}>{mut1}, {pos2}{ref2}>{mut2})"
                        double_mutations[double_key] += count
        
        if not mutation_types:
            return ""
            
        # Create single mutation spectrum
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        types = sorted(mutation_types.keys())
        if not types:
            return ""
            
        counts = [mutation_types[t] for t in types]
        if not counts or total_mutations == 0:
            return ""
            
        percentages = [100 * c / total_mutations for c in counts]
        
        sns.barplot(x=types, y=percentages)
        plt.title('Single Mutation Spectrum')
        plt.xlabel('Mutation Type')
        plt.ylabel('Percentage of Total Mutations')
        plt.xticks(rotation=45)
        
        # Create linked mutation spectrum (top 10 most frequent)
        plt.subplot(1, 2, 2)
        if double_mutations:
            sorted_doubles = sorted(double_mutations.items(), key=lambda x: x[1], reverse=True)[:10]
            if sorted_doubles:
                double_types, double_counts = zip(*sorted_doubles)
                double_percentages = [100 * c / total_mutations for c in double_counts]
                
                sns.barplot(x=list(range(len(double_types))), y=double_percentages)
                plt.title('Top 10 Linked Mutations')
                plt.xlabel('Mutation Pair')
                plt.ylabel('Percentage of Total Mutations')
                plt.xticks(range(len(double_types)), double_types, rotation=45, ha='right')
        
        plt.tight_layout()
        return fig_to_base64()
    except Exception as e:
        logger.error(f"Error creating mutation spectrum: {str(e)}")
        return ""

def create_position_mutation_plot(results_df: pd.DataFrame, reference_seq: str) -> str:
    """Create a plot showing number of different mutations at each position."""
    if results_df.empty:
        return ""
    
    # Track unique mutations at each position
    position_mutations = defaultdict(set)
    total_reads = results_df['count'].sum()
    
    if total_reads == 0:
        return ""
    
    # Process each haplotype
    for _, row in results_df.iterrows():
        haplotype = row['haplotype']
        for pos, (ref, var) in enumerate(zip(reference_seq.upper(), haplotype)):
            if var.islower():  # This is a mutation
                mutation = f"{ref}>{var.upper()}"
                position_mutations[pos + 1].add(mutation)  # 1-based position
    
    if not position_mutations:
        return ""
    
    # Convert to plot data
    positions = sorted(position_mutations.keys())
    if not positions:
        return ""
        
    mutation_counts = [len(position_mutations[pos]) for pos in positions]
    if not mutation_counts:
        return ""
    
    # Create Plotly figure
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=positions,
        y=mutation_counts,
        marker_color='rgb(158,202,225)',
        hovertemplate="Position: %{x}<br>" +
                     "Different mutations: %{y}<br>" +
                     "Mutations: %{customdata}<extra></extra>",
        customdata=[list(position_mutations[pos]) for pos in positions]
    ))
    
    max_count = max(mutation_counts) if mutation_counts else 0
    
    fig.update_layout(
        title='Mutation Diversity by Position',
        xaxis=dict(
            title='Position in Sequence',
            tickmode='linear'
        ),
        yaxis=dict(
            title='Number of Different Mutations',
            range=[0, max_count * 1.1] if max_count > 0 else [0, 1]
        ),
        margin=dict(b=50, l=60, r=20, t=40),
        width=800,
        height=400,
        template='plotly_white'
    )
    
    # Convert to HTML
    plot_html = fig.to_html(
        full_html=False,
        include_plotlyjs='cdn',
        config={'displayModeBar': True}
    )
    
    return plot_html

def create_indel_summary(results_df: pd.DataFrame) -> pd.DataFrame:
    """Create a summary of indels found in the results.
    
    Args:
        results_df: DataFrame containing mutation results
        
    Returns:
        DataFrame with indel summary statistics
    """
    if results_df.empty:
        return pd.DataFrame()
        
    indel_data = []
    for _, row in results_df.iterrows():
        haplotype = row['haplotype']
        count = row['count']
        frequency = row['frequency']
        
        # Look for indels (marked by '-' in the haplotype)
        for pos, base in enumerate(haplotype):
            if base == '-':
                indel_data.append({
                    'position': pos + 1,
                    'type': 'deletion',
                    'size': 1,
                    'ref': 'N',  # We don't have the reference base in this context
                    'alt': '-',
                    'in_bed': False,  # We don't have BED information in this context
                    'frequency': frequency,
                    'count': count
                })
        for pos, bases in AmpliconProcessor.parse_insertions(row.get('insertions', '')):
            indel_data.append({
                'position': pos,
                'type': 'insertion',
                'size': len(bases),
                'ref': '-',
                'alt': bases,
                'in_bed': False,
                'frequency': frequency,
                'count': count
            })
    
    if not indel_data:
        return pd.DataFrame()
        
    indel_df = pd.DataFrame(indel_data)
    
    # Create summary table with rich formatting
    summary = []
    for (type_name, in_bed), group in indel_df.groupby(['type', 'in_bed']):
        summary.append({
            'Type': type_name.capitalize(),
            'In BED': 'Yes' if in_bed else 'No',
            'Count': len(group),
            'Avg Size': f"{group['size'].mean():.1f}",
            'Max Size': group['size'].max(),
            'Avg Frequency': f"{group['frequency'].mean():.2%}",
            'Positions': ', '.join(map(str, sorted(group['position'].unique())))
        })
    
    return pd.DataFrame(summary)

def create_indel_plot(results_df: pd.DataFrame) -> str:
    """Create indel distribution plot."""
    if results_df.empty:
        return ""
        
    # Extract indel information
    indel_data = []
    total_reads = results_df['count'].sum()
    
    if total_reads == 0:
        return ""
    
    for _, row in results_df.iterrows():
        haplotype = row['haplotype']
        count = row['count']
        
        # Find indels (marked by '-')
        current_indel = []
        for i, base in enumerate(haplotype):
            if base == '-':
                current_indel.append(i)
            elif current_indel:
                indel_data.append({
                    'position': min(current_indel),
                    'size': len(current_indel),
                    'count': count
                })
                current_indel = []
                
        # Handle indel at end of sequence
        if current_indel:
            indel_data.append({
                'position': min(current_indel),
                'size': len(current_indel),
                'count': count
            })

        for pos, bases in AmpliconProcessor.parse_insertions(row.get('insertions', '')):
            indel_data.append({
                'position': pos,
                'size': len(bases),
                'count': count
            })
    
    if not indel_data:
        return ""
    
    # Convert to DataFrame for easier plotting
    indel_df = pd.DataFrame(indel_data)
    
    # Calculate frequencies
    indel_df['frequency'] = indel_df['count'] / total_reads * 100
    
    # Create size distribution plot
    fig = go.Figure()
    
    # Group by size and sum frequencies
    size_dist = indel_df.groupby('size')['frequency'].sum().reset_index()
    if size_dist.empty:
        return ""
    
    fig.add_trace(go.Bar(
        x=size_dist['size'],
        y=size_dist['frequency'],
        name='Indel Size Distribution',
        marker_color='rgb(55, 83, 109)'
    ))
    
    fig.update_layout(
        title='Indel Size Distribution',
        xaxis=dict(
            title='Indel Size (bp)',
            tickmode='linear'
        ),
        yaxis=dict(
            title='Frequency (%)',
            range=[0, size_dist['frequency'].max() * 1.1] if not size_dist.empty else [0, 100]
        ),
        showlegend=False,
        margin=dict(l=60, r=20, t=40, b=60),
        width=800,
        height=400,
        template='plotly_white'
    )
    
    # Convert to HTML
    plot_html = fig.to_html(
        full_html=False,
        include_plotlyjs='cdn',
        config={'displayModeBar': True}
    )
    
    return plot_html

def create_rv0678_analysis_table(results_df: pd.DataFrame) -> pd.DataFrame:
    """Create a detailed analysis table for rv0678 data."""
    if results_df.empty:
        return pd.DataFrame()
    
    # Calculate mutation statistics
    total_reads = results_df['count'].sum()
    unique_haplotypes = len(results_df)
    
    # Count unique single mutations more accurately
    single_mutation_haplotypes = results_df[results_df['mutations'] == 1]
    unique_single_mutations = len(single_mutation_haplotypes)
    
    # Calculate theoretical maximum single mutations
    reference_length = len(single_mutation_haplotypes.iloc[0]['haplotype']) if not single_mutation_haplotypes.empty else 0
    theoretical_max_single_mutations = reference_length * 3  # 3 possible mutations per position (A->C,G,T, etc.)
    
    max_frequency = (results_df['count'].max() / total_reads * 100)
    avg_mutations = (results_df['mutations'].astype(float) * results_df['count']).sum() / total_reads
    full_length_reads = results_df[results_df['is_full_length']]['count'].sum()
    full_length_percent = (full_length_reads / total_reads * 100)
    
    # Create mutation distribution
    mutation_dist = results_df.groupby('mutations')['count'].sum()
    mutation_dist_percent = (mutation_dist / total_reads * 100).round(2)
    
    # Create analysis table
    analysis_data = {
        'Metric': [
            'Total Reads',
            'Unique Haplotypes',
            'Unique Single Mutations',
            'Theoretical Max Single Mutations',
            'Maximum Frequency (%)',
            'Average Mutations',
            'Full Length Reads',
            'Full Length Percentage (%)',
            'Number of References'
        ],
        'Value': [
            f"{total_reads:,}",
            f"{unique_haplotypes:,}",
            f"{unique_single_mutations:,}",
            f"{theoretical_max_single_mutations:,}",
            f"{max_frequency:.2f}",
            f"{avg_mutations:.2f}",
            f"{full_length_reads:,}",
            f"{full_length_percent:.2f}",
            "1"
        ]
    }
    
    return pd.DataFrame(analysis_data)

def format_summary_table(results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Format complete summary table with all mutation statistics."""
    summary_data = []
    
    for sample, df in results.items():
        if df.empty:
            continue
            
        # Basic stats
        total_reads = df['count'].sum()
        unique_haplotypes = len(df)
        max_freq = (df['count'].max() / total_reads * 100) if total_reads > 0 else 0.0
        avg_mutations = (df['mutations'].astype(float) * df['count']).sum() / total_reads if total_reads > 0 else 0.0
        
        # Full length statistics
        full_length_reads = df[df['is_full_length']]['count'].sum()
        full_length_percent = (full_length_reads / total_reads * 100) if total_reads > 0 else 0.0
        
        # Single mutation statistics
        single_mut_reads = df[df['mutations'] == 1]['count'].sum()
        single_mut_percent = (single_mut_reads / total_reads * 100) if total_reads > 0 else 0.0
        
        # Full length single mutations
        full_length_single = df[(df['mutations'] == 1) & (df['is_full_length'])]['count'].sum()
        full_length_single_percent = (full_length_single / total_reads * 100) if total_reads > 0 else 0.0
        
        summary_data.append({
            'sample': sample,
            'total_reads': total_reads,
            'unique_haplotypes': unique_haplotypes,
            'max_frequency': max_freq,
            'avg_mutations': avg_mutations,
            'full_length_reads': full_length_reads,
            'full_length_percent': full_length_percent,
            'single_mutations': single_mut_reads,
            'single_mutation_percent': single_mut_percent,
            'full_length_single_mutations': full_length_single,
            'full_length_single_percent': full_length_single_percent
        })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Format for display
    if not summary_df.empty:
        # Round percentages to 2 decimal places
        percentage_cols = ['max_frequency', 'full_length_percent', 
                         'single_mutation_percent', 'full_length_single_percent']
        summary_df[percentage_cols] = summary_df[percentage_cols].round(2)
        
        # Round avg_mutations to 2 decimal places
        summary_df['avg_mutations'] = summary_df['avg_mutations'].round(2)
    
    return summary_df


def generate_report(results: Dict[str, pd.DataFrame], 
                   summary: pd.DataFrame, 
                   output_path: Path,
                   reference_seq: str):
    """Generate HTML report with analysis results."""
    try:
        # Ensure we have results to process
        if not results or all(df.empty for df in results.values()):
            logger.warning("No results to generate report from")
            with open(output_path, 'w') as f:
                f.write("""
                <html>
                <head><title>CloneArmy Analysis Report - No Results</title></head>
                <body>
                <h1>CloneArmy Analysis Report</h1>
                <p>No valid results were found to generate a report. This could be due to:</p>
                <ul>
                    <li>No reads meeting quality thresholds</li>
                    <li>No haplotypes meeting minimum read count</li>
                    <li>All reads filtered out during processing</li>
                </ul>
                </body>
                </html>
                """)
            return

        # Create plots
        mutation_freq_plot = ""
        mutation_spectrum = ""
        position_plot = ""
        indel_plot = ""
        
        # Only create plots if we have valid data
        for ref_name, df in results.items():
            if not df.empty and len(df) > 0:
                mutation_freq_plot = create_mutation_frequency_plot(df, reference_seq)
                mutation_spectrum = create_mutation_spectrum({ref_name: df}, reference_seq)
                position_plot = create_position_mutation_plot(df, reference_seq)
                indel_plot = create_indel_plot(df)
                break

        # Format summary table
        summary_table = format_summary_table(results)
        
        # Get template path
        template_path = Path(__file__).parent / 'templates' / 'report_template.html'
        if not template_path.exists():
            # Inline template if file doesn't exist
            template_str = """
            <!DOCTYPE html>
            <html>
            <head>
                <title>CloneArmy Analysis Report</title>
                <style>
                    body { font-family: Arial, sans-serif; margin: 20px; }
                    table { border-collapse: collapse; margin: 10px 0; }
                    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                    th { background-color: #f5f5f5; }
                    .plot { margin: 20px 0; }
                    .warning { color: #856404; background-color: #fff3cd; padding: 10px; border-radius: 4px; }
                </style>
            </head>
            <body>
                <h1>CloneArmy Analysis Report</h1>
                <p>Generated on: {{ timestamp }}</p>
                
                {% if summary_table is not none %}
                <h2>Summary</h2>
                {{ summary_table | safe }}
                {% endif %}
                
                {% if mutation_freq_plot %}
                <h2>Mutation Frequency</h2>
                <div class="plot">{{ mutation_freq_plot | safe }}</div>
                {% endif %}
                
                {% if mutation_spectrum %}
                <h2>Mutation Spectrum</h2>
                <div class="plot">
                    <img src="{{ mutation_spectrum }}" alt="Mutation Spectrum">
                </div>
                {% endif %}
                
                {% if position_plot %}
                <h2>Mutations by Position</h2>
                <div class="plot">{{ position_plot | safe }}</div>
                {% endif %}
                
                {% if indel_plot %}
                <h2>Indel Distribution</h2>
                <div class="plot">{{ indel_plot | safe }}</div>
                {% endif %}
            </body>
            </html>
            """
        else:
            with open(template_path) as f:
                template_str = f.read()
        
        # Render template
        template = Template(template_str)
        html = template.render(
            timestamp=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            summary_table=summary_table.to_html(index=False) if not summary_table.empty else None,
            mutation_freq_plot=mutation_freq_plot,
            mutation_spectrum=mutation_spectrum,
            position_plot=position_plot,
            indel_plot=indel_plot
        )
        
        # Write report
        with open(output_path, 'w') as f:
            f.write(html)
            
    except Exception as e:
        logger.error(f"Error generating report: {str(e)}")
        raise

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Mutation Analysis Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body { padding: 20px; }
        .plot-container { margin-bottom: 30px; }
        .plot-description { 
            font-style: italic;
            color: #666;
            margin-top: 10px;
        }
        .stats-box {
            background-color: #f8f9fa;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
        }
        .table-responsive {
            margin-bottom: 30px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1 class="mb-4">Mutation Analysis Report</h1>
        <p class="text-muted">Generated on {{ date }}</p>

        {% if has_data %}
            <!-- rv0678 Analysis Section -->
            <div class="card mb-4">
                <div class="card-header">
                    <h2>rv0678 Gene Analysis</h2>
                </div>
                <div class="card-body">
                    <div class="table-responsive">
                        {{ rv0678_analysis | safe }}
                    </div>
                </div>
            </div>

            <!-- Summary Section -->
            <div class="card mb-4">
                <div class="card-header">
                    <h2>Summary</h2>
                </div>
                <div class="card-body">
                    <div class="table-responsive">
                        {{ summary | safe }}
                    </div>
                </div>
            </div>

            <!-- Indel Analysis Section -->
            <div class="card mb-4">
                <div class="card-header">
                    <h2>Indel Analysis</h2>
                </div>
                <div class="card-body">
                    <h3>Indel Summary</h3>
                    <div class="table-responsive">
                        {{ indel_summary | safe }}
                    </div>
                    
                    {% if indel_plot %}
                    <h3 class="mt-4">Indel Distribution</h3>
                    <div class="plot-container">
                        {{ indel_plot | safe }}
                        <p class="plot-description">
                            Distribution of insertions and deletions across the sequence.
                            Stars indicate indels that overlap with regions in the provided BED file.
                            Hover over points to see detailed information.
                        </p>
                    </div>
                    {% endif %}
                </div>
            </div>

            <!-- Mutation Analysis Section -->
            <div class="card mb-4">
                <div class="card-header">
                    <h2>Mutation Analysis</h2>
                </div>
                <div class="card-body">
                    <!-- Mutation Spectrum -->
                    {% if mutation_spectrum %}
                    <div class="plot-container">
                        {{ mutation_spectrum | safe }}
                        <p class="plot-description">
                            Distribution of mutation types across all samples.
                            Shows the frequency of different types of base changes.
                        </p>
                    </div>
                    {% endif %}
                    
                    <!-- Mutation Frequency Plot -->
                    {% if mutation_freq_plot %}
                    <div class="plot-container">
                        {{ mutation_freq_plot | safe }}
                        <p class="plot-description">
                            Distribution of individual mutation frequencies across the sequence.
                            Each point represents a unique mutation and its frequency.
                        </p>
                    </div>
                    {% endif %}
                    
                    <!-- Position Mutation Plot -->
                    {% if position_plot %}
                    <div class="plot-container">
                        {{ position_plot | safe }}
                        <p class="plot-description">
                            Number of different mutations observed at each position.
                            Hover over bars to see the specific mutations at each position.
                        </p>
                    </div>
                    {% endif %}
                </div>
            </div>

            <!-- Statistics Section -->
            <div class="card mb-4">
                <div class="card-header">
                    <h2>Statistics</h2>
                </div>
                <div class="card-body">
                    {% for sample, stat in stats.items() %}
                    <div class="stats-box">
                        <h3>{{ sample }}</h3>
                        <ul class="list-unstyled">
                            <li><strong>Total Reads:</strong> {{ "{:,}".format(stat.total_reads) }}</li>
                            <li><strong>Unique Haplotypes:</strong> {{ "{:,}".format(stat.unique_haplotypes) }}</li>
                            <li><strong>Full Length Reads:</strong> {{ "{:.2f}%".format(stat.full_length_percent) }}</li>
                        </ul>
                        <h4>Mutation Rates:</h4>
                        <ul class="list-unstyled">
                        {% for mutations, rate in stat.mutation_rates %}
                            <li>{{ mutations }} mutations: {{ "{:.2f}%".format(rate) }}</li>
                        {% endfor %}
                        </ul>
                    </div>
                    {% endfor %}
                </div>
            </div>
        {% else %}
            <div class="alert alert-warning">
                No mutation data available for analysis.
            </div>
        {% endif %}
    </div>
</body>
</html>
"""