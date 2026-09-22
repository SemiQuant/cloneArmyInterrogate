from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import numpy as np
from typing import Dict, Tuple, List, Set, Union
import logging
from pathlib import Path
import plotly.graph_objects as go
from collections import defaultdict

from .processor import AmpliconProcessor

logger = logging.getLogger(__name__)

def get_mutations(haplotype: str, reference: str, insertions: str = '') -> Set[Tuple[int, str, str]]:
    """Get all mutations (including insertions) in a haplotype compared to reference."""
    mutations = set()
    positions = []
    ref_bases = []
    mut_bases = []
    
    # First collect all mutations
    for pos, (ref_base, hap_base) in enumerate(zip(reference, haplotype)):
        if (hap_base.islower() or      # Regular mutation (lowercase)
            hap_base == '-' or         # Deletion
            ref_base == '-'):          # Insertion
            positions.append(pos + 1)
            ref_bases.append(ref_base.upper())
            mut_bases.append(hap_base.upper())

    for pos, bases in AmpliconProcessor.parse_insertions(insertions):
        positions.append(pos)
        ref_bases.append('-')
        mut_bases.append(f"+{bases}")
    
    # Add single mutations
    for i in range(len(positions)):
        mutations.add((positions[i], ref_bases[i], mut_bases[i]))
    
    # Add linked double mutations
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            mutations.add((
                positions[i], positions[j],  # positions
                ref_bases[i], ref_bases[j],  # reference bases
                mut_bases[i], mut_bases[j]   # mutant bases
            ))
    
    return mutations

def perform_proportion_test(count1: int, total1: int, count2: int, total2: int) -> float:
    """Perform chi-square test of independence for two proportions."""
    contingency = np.array([
        [count1, total1 - count1],
        [count2, total2 - count2]
    ])
    _, pval, _, _ = stats.chi2_contingency(contingency)
    return pval

def compare_samples(results1: pd.DataFrame, 
                   results2: pd.DataFrame,
                   ref_seq: str,
                   full_length_only: bool = False) -> pd.DataFrame:
    """
    Compare mutation and haplotype frequencies between two samples.
    
    Args:
        results1: Results DataFrame for first sample
        results2: Results DataFrame for second sample
        ref_seq: Reference sequence string
        full_length_only: Only consider sequences that cover the entire reference
    """
    # Filter for full-length sequences if requested
    if full_length_only:
        results1 = results1[results1['is_full_length'] == True].copy()
        results2 = results2[results2['is_full_length'] == True].copy()
        
        if results1.empty or results2.empty:
            logger.warning("No full-length sequences found in one or both samples")
            return pd.DataFrame()
    
    # Get total reads for each sample
    total_reads1 = results1['count'].sum()
    total_reads2 = results2['count'].sum()
    
    logger.info(f"Sample 1: {total_reads1:,} reads" + 
                " (full-length only)" if full_length_only else "")
    logger.info(f"Sample 2: {total_reads2:,} reads" + 
                " (full-length only)" if full_length_only else "")
    
    # Process each haplotype and its mutations
    mutations1 = defaultdict(int)
    mutations2 = defaultdict(int)
    
    # Get mutations for each sample
    for df, mutations_dict in [(results1, mutations1), (results2, mutations2)]:
        for _, row in df.iterrows():
            haplotype = row['haplotype']
            count = row['count']
            for mut in get_mutations(haplotype, ref_seq, row.get('insertions', '')):
                mutations_dict[mut] += count
    
    # Compare frequencies for each mutation
    all_mutations = set(mutations1.keys()) | set(mutations2.keys())
    comparisons = []
    
    for mut in all_mutations:
        if len(mut) == 3:  # Single mutation
            pos, ref_base, mut_base = mut
            mutation_type = 'substitution'
            position_str = str(pos)
            ref_str = ref_base
            mut_str = mut_base
        else:  # Linked double mutation
            pos1, pos2, ref1, ref2, mut1, mut2 = mut
            mutation_type = 'double'
            position_str = f"{pos1},{pos2}"
            ref_str = f"{ref1},{ref2}"
            mut_str = f"{mut1},{mut2}"
        
        count1 = mutations1.get(mut, 0)
        count2 = mutations2.get(mut, 0)
        
        # Calculate frequencies
        freq1 = (count1 / total_reads1) * 100 if total_reads1 > 0 else 0
        freq2 = (count2 / total_reads2) * 100 if total_reads2 > 0 else 0
        
        # Perform statistical test
        pval = perform_proportion_test(count1, total_reads1, count2, total_reads2)
        
        # Store comparison results
        comparisons.append({
            'reference_name': results1['reference'].iloc[0],
            'position': position_str,
            'reference_base': ref_str,
            'mutation': mut_str,
            'mutation_type': mutation_type,
            'sample1_percent': freq1,
            'sample1_reads': count1,
            'sample2_percent': freq2,
            'sample2_reads': count2,
            'pvalue': pval
        })
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(comparisons)
    if not results_df.empty:
        results_df['fdr_pvalue'] = fdrcorrection(results_df['pvalue'])[1]
        results_df = results_df.sort_values(['position', 'fdr_pvalue'])
    
    return results_df

def create_mutation_plot(comparison_df: pd.DataFrame, output_path: Path) -> None:
    """Create a Plotly plot showing mutation frequencies."""
    # Sort by numeric position
    comparison_df['position_num'] = comparison_df['position'].astype(str).str.extract('(\d+)').astype(int)
    comparison_df = comparison_df.sort_values('position_num')
    
    # Find the percentage columns - they should end with '_percent'
    percent_cols = [col for col in comparison_df.columns if col.endswith('_percent')]
    if len(percent_cols) < 2:
        logger.error(f"Expected two percentage columns, found: {percent_cols}")
        return None
    
    # Create x-axis labels
    comparison_df['x_label'] = comparison_df.apply(
        lambda x: f"{x['position']} {x['mutation']}", 
        axis=1
    )
    
    # Create the figure
    fig = go.Figure()
    
    # Add first sample percentages
    fig.add_trace(go.Bar(
        x=comparison_df['x_label'],
        y=comparison_df[percent_cols[0]],
        name=percent_cols[0].replace('_percent', ''),
        marker_color='rgb(31, 119, 180)'
    ))
    
    # Add second sample percentages (negative values)
    fig.add_trace(go.Bar(
        x=comparison_df['x_label'],
        y=-comparison_df[percent_cols[1]],
        name=percent_cols[1].replace('_percent', ''),
        marker_color='rgb(255, 127, 14)'
    ))
    
    # Update layout
    fig.update_layout(
        barmode='relative',
        template='plotly_white',
        xaxis=dict(
            title='Position and Mutation',
            tickangle=45,
            tickfont=dict(size=10),
            categoryorder='array',  # Force the x-axis to use our sorted order
            categoryarray=comparison_df['x_label']
        ),
        yaxis=dict(
            title='Read Percent',
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='black',
        ),
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='right',
            x=1
        ),
        margin=dict(
            l=50,
            r=50,
            t=50,
            b=100
        )
    )
    
    # Save plot
    fig.write_html(output_path)
    logger.info(f"Plot saved to: {output_path}")

def run_comparative_analysis(
    results1: Dict[str, pd.DataFrame],
    results2: Dict[str, pd.DataFrame],
    output_dir: Union[str, Path],
    full_length_only: bool = False
) -> pd.DataFrame:
    """
    Run comparative analysis between two sets of results.
    
    Args:
        results1: Results from first set of samples
        results2: Results from second set of samples
        output_dir: Directory for output files
        full_length_only: Only consider sequences that cover the entire reference
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get common references between the two sets
    refs1 = set()
    refs2 = set()
    for df in results1.values():
        refs1.update(df['reference'].unique())
    for df in results2.values():
        refs2.update(df['reference'].unique())
    common_refs = refs1 & refs2
    
    if not common_refs:
        logger.warning("No common references found between the two sets")
        return pd.DataFrame()
    
    # Process each reference
    all_comparisons = []
    for ref in common_refs:
        # Get results for this reference from each set
        ref_results1 = pd.concat([
            df[df['reference'] == ref] for df in results1.values()
        ])
        ref_results2 = pd.concat([
            df[df['reference'] == ref] for df in results2.values()
        ])
        
        if ref_results1.empty or ref_results2.empty:
            continue
            
        # Get reference sequence from first processor
        ref_seq = ref_results1['reference'].iloc[0]
        
        # Compare samples
        comparison_df = compare_samples(
            ref_results1, ref_results2,
            ref_seq=ref_seq,
            full_length_only=full_length_only
        )
        
        if not comparison_df.empty:
            all_comparisons.append(comparison_df)
    
    if not all_comparisons:
        return pd.DataFrame()
    
    # Combine all comparisons
    final_df = pd.concat(all_comparisons, ignore_index=True)
    
    # Create plots
    for ref in common_refs:
        ref_df = final_df[final_df['reference_name'] == ref]
        if not ref_df.empty:
            plot_path = output_dir / f'mutation_plot_{ref}.html'
            create_mutation_plot(ref_df, plot_path)
    
    # Save results
    final_df.to_csv(output_dir / 'comparison_results.csv', index=False)
    
    return final_df