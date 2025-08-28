import marimo

__generated_with = "0.13.15"
app = marimo.App(
    width="medium",
    layout_file="layouts/AIM_Analysis_for_Anopheles_Species.grid.json",
)


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Ancestry analysis using ancestry informative markers (AIMs) and Inversions karyotyping 

    In this notebook we will load a vcf and plot the ancestry components of samples using ancestry markers
    """
    )
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import zarr
    import allel
    from collections import defaultdict
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns
    import marimo as mo
    from pathlib import Path

    print("All dependencies imported successfully!")
    print("Available modules:")
    print(f"- pandas: {pd.__version__}")
    print(f"- numpy: {np.__version__}")
    print(f"- scikit-allel: {allel.__version__}")
    return allel, mo, np, os, pd, plt, sns, zarr


@app.cell
def _(os):
    # Your file paths
    multisample_vcf_path = "/Users/talal.alyazeedi/Documents/LSTM_Postdoc/A.gambiae_BSA/pooled.continuous.freebayes.parallel.vcf.gz"
    aims_zarr_gambcolu_path = "/Users/talal.alyazeedi/Documents/LSTM_Postdoc/A.gambiae_BSA/rna-seq-pop/resources/gamb_vs_colu.zarr"
    aims_zarr_arab_path = "/Users/talal.alyazeedi/Documents/LSTM_Postdoc/A.gambiae_BSA/rna-seq-pop/resources/gambcolu_vs_arab.zarr"

    # Output directory
    output_directory = "/Users/talal.alyazeedi/Documents/LSTM_Postdoc/A.gambiae_BSA/ancestry_analysis_results"

    # Analysis parameters
    target_contigs = ['2L', '2R', '3L', '3R', 'X']  # Standard Anopheles chromosomes
    quality_filter = 30          # Quality filter threshold
    missing_proportion = 0.1     # Missing data proportion (10%)
    organism_ploidy = 2         # Diploid

    # Create output directory
    os.makedirs(output_directory, exist_ok=True)

    print("Configuration loaded with your paths:")
    print(f"VCF file: {multisample_vcf_path}")
    print(f"AIMs (gambiae/coluzzii): {aims_zarr_gambcolu_path}")
    print(f"AIMs (arabiensis): {aims_zarr_arab_path}")
    print(f"Output directory: {output_directory}")
    print(f"Chromosomes to analyze: {target_contigs}")
    return (
        aims_zarr_arab_path,
        aims_zarr_gambcolu_path,
        multisample_vcf_path,
        target_contigs,
    )


@app.cell
def _(allel, multisample_vcf_path, np, target_contigs):
    print("Loading VCF file...")

    vcf_dataset = allel.read_vcf(multisample_vcf_path)

    if vcf_dataset is not None:
        print("VCF loaded successfully!")
        print("VCF data structure:")
        for data_key in sorted(vcf_dataset.keys()):
            if hasattr(vcf_dataset[data_key], 'shape'):
                print(f"  {data_key}: shape {vcf_dataset[data_key].shape}")
            else:
                print(f"  {data_key}: {type(vcf_dataset[data_key])}")

        # Get basic info
        if 'samples' in vcf_dataset:
            vcf_sample_names = vcf_dataset['samples']
            print(f"\nSamples ({len(vcf_sample_names)}): {list(vcf_sample_names)}")
        if 'variants/POS' in vcf_dataset:
            vcf_positions = vcf_dataset['variants/POS']
            print(f"\nTotal variants: {len(vcf_positions):,}")
        if 'variants/CHROM' in vcf_dataset:
            vcf_chromosomes = vcf_dataset['variants/CHROM']
            unique_chromosomes = np.unique(vcf_chromosomes)
            print(f"\nChromosomes found: {unique_chromosomes}")

            # Count variants per chromosome
            print("Variants per chromosome:")
            for chromosome in target_contigs:
                variant_count = np.sum(vcf_chromosomes == chromosome)
                print(f"  {chromosome}: {variant_count:,} variants")
        # Check genotype data
        if 'calldata/GT' in vcf_dataset:
            vcf_genotypes = vcf_dataset['calldata/GT']
            print(f"\nGenotypes shape: {vcf_genotypes.shape}")
            print(f"   (variants × samples × ploidy)")
    else:
        print("Failed to load VCF file")
    return (
        vcf_chromosomes,
        vcf_dataset,
        vcf_genotypes,
        vcf_positions,
        vcf_sample_names,
    )


@app.cell
def _(aims_zarr_arab_path, aims_zarr_gambcolu_path, target_contigs, zarr):
    print("Loading AIMs data...")

    # Load gambiae vs coluzzii AIMs
    aims_gambcolu = zarr.open(aims_zarr_gambcolu_path, mode='r')
    print("Gambiae/Coluzzii AIMs data loaded successfully!")
    print(f"Available chromosomes: {list(aims_gambcolu.keys())}")

    # Check each chromosome
    for contig in target_contigs:
        if contig in aims_gambcolu:
            n_aims = len(aims_gambcolu[contig]['POS'])
            print(f"   {contig}: {n_aims} AIMs")
        else:
            print(f"   {contig}: Not found")

    # Load arabiensis AIMs
    try:
        aims_arab = zarr.open(aims_zarr_arab_path, mode='r')
        print(f"\nArabiensis AIMs data loaded successfully!")
        print(f"Available chromosomes: {list(aims_arab.keys())}")

        # Check each chromosome for arabiensis
        for contig_name in target_contigs:
            if contig_name in aims_arab:
                n_aims_arab = len(aims_arab[contig_name]['POS'])
                print(f"   {contig_name}: {n_aims_arab} arabiensis AIMs")
            else:
                print(f"   {contig_name}: Not found")

    except Exception as e:
        print(f"Could not load arabiensis AIMs: {e}")
        aims_arab = None
    return aims_arab, aims_gambcolu


@app.cell
def _(
    aims_arab,
    aims_gambcolu,
    target_contigs,
    vcf_chromosomes,
    vcf_positions,
):
    print("Finding AIMs overlap with VCF data...")

    # Initialize results dictionaries
    gambcolu_overlap_results = {}
    arab_overlap_results = {}

    # Check gambiae vs coluzzii AIMs overlap
    print("\n=== Gambiae vs Coluzzii AIMs Analysis ===")
    for chrom_id in target_contigs:
        print(f"\nAnalyzing chromosome {chrom_id}...")

        # Check if chromosome exists in VCF
        chrom_mask = vcf_chromosomes == chrom_id
        chrom_vcf_positions = vcf_positions[chrom_mask]

        if len(chrom_vcf_positions) == 0:
            print(f"  No variants found for {chrom_id} in VCF")
            continue

        print(f"  VCF: {len(chrom_vcf_positions):,} variants")

        # Get AIMs data
        if chrom_id in aims_gambcolu:
            chrom_aims_positions = aims_gambcolu[chrom_id]['POS'][:]
            print(f"  AIMs: {len(chrom_aims_positions)} markers")

            # Find intersection
            vcf_pos_set = set(chrom_vcf_positions)
            aims_pos_set = set(chrom_aims_positions)
            common_positions = vcf_pos_set.intersection(aims_pos_set)

            overlap_percentage = len(common_positions)/len(chrom_aims_positions)*100
            print(f"  Overlap: {len(common_positions)} positions ({overlap_percentage:.1f}% of AIMs)")

            if len(common_positions) > 0:
                gambcolu_overlap_results[chrom_id] = {
                    'vcf_positions': chrom_vcf_positions,
                    'aims_positions': chrom_aims_positions,
                    'common_positions': common_positions,
                    'overlap_count': len(common_positions)
                }
        else:
            print(f"  No AIMs data for {chrom_id}")

    # Check arabiensis AIMs overlap if available
    if aims_arab is not None:
        print("\n=== Arabiensis vs Gambiae/Coluzzii AIMs Analysis ===")
        for chr_id in target_contigs:
            print(f"\nAnalyzing chromosome {chr_id}...")

            # Get VCF positions for this chromosome
            chr_mask = vcf_chromosomes == chr_id
            chr_vcf_positions = vcf_positions[chr_mask]

            if len(chr_vcf_positions) == 0:
                continue

            print(f"  VCF: {len(chr_vcf_positions):,} variants")

            # Get arabiensis AIMs data
            if chr_id in aims_arab:
                chr_aims_positions = aims_arab[chr_id]['POS'][:]
                print(f"  Arabiensis AIMs: {len(chr_aims_positions)} markers")

                # Find intersection
                vcf_positions_set = set(chr_vcf_positions)
                aims_positions_set = set(chr_aims_positions)
                overlapping_positions = vcf_positions_set.intersection(aims_positions_set)

                overlap_pct = len(overlapping_positions)/len(chr_aims_positions)*100
                print(f"  Overlap: {len(overlapping_positions)} positions ({overlap_pct:.1f}% of AIMs)")

                if len(overlapping_positions) > 0:
                    arab_overlap_results[chr_id] = {
                        'vcf_positions': chr_vcf_positions,
                        'aims_positions': chr_aims_positions,
                        'common_positions': overlapping_positions,
                        'overlap_count': len(overlapping_positions)
                    }
            else:
                print(f"  No arabiensis AIMs data for {chr_id}")

    # Summary
    print(f"\n=== Summary ===")
    print(f"Gambiae/Coluzzii - Chromosomes with overlap: {list(gambcolu_overlap_results.keys())}")
    total_gambcolu_aims = sum([result['overlap_count'] for result in gambcolu_overlap_results.values()])
    print(f"Total gambiae/coluzzii overlapping AIMs: {total_gambcolu_aims}")

    if arab_overlap_results:
        print(f"Arabiensis - Chromosomes with overlap: {list(arab_overlap_results.keys())}")
        total_arab_aims = sum([result['overlap_count'] for result in arab_overlap_results.values()])
        print(f"Total arabiensis overlapping AIMs: {total_arab_aims}")
    else:
        print("No arabiensis overlap results")
    return arab_overlap_results, gambcolu_overlap_results


@app.cell
def _(
    aims_arab,
    aims_gambcolu,
    arab_overlap_results,
    gambcolu_overlap_results,
    np,
    vcf_chromosomes,
    vcf_dataset,
    vcf_genotypes,
    vcf_positions,
):
    def analyze_ancestry_for_dataset(overlap_results, aims_dataset, dataset_name):
        """Analyze ancestry for a given dataset (gambcolu or arab)"""

        print(f"\nAnalyzing {dataset_name} ancestry...")
        ancestry_results = {}

        for chr_name in overlap_results.keys():
            print(f"Processing chromosome {chr_name}...")

            # Get VCF data for this chromosome
            chr_mask = vcf_chromosomes == chr_name
            chr_vcf_pos = vcf_positions[chr_mask]
            chr_vcf_geno = vcf_genotypes[chr_mask]
            chr_vcf_ref = vcf_dataset['variants/REF'][chr_mask]
            chr_vcf_alt = vcf_dataset['variants/ALT'][chr_mask]

            # Get AIMs data
            chr_aims_pos = aims_dataset[chr_name]['POS'][:]
            if dataset_name == "gambcolu":
                chr_aims_anc1 = aims_dataset[chr_name]['gamb_allele'][:]  # gambiae
                chr_aims_anc2 = aims_dataset[chr_name]['colu_allele'][:]  # coluzzii
            else:  # arabiensis
                chr_aims_anc1 = aims_dataset[chr_name]['gambcolu_allele'][:]  # gambiae+coluzzii
                chr_aims_anc2 = aims_dataset[chr_name]['arab_allele'][:]     # arabiensis

            # Find matching positions
            common_pos = overlap_results[chr_name]['common_positions']

            # Get indices for matching positions
            vcf_indices = []
            aims_indices = []

            for pos in common_pos:
                vcf_idx = np.where(chr_vcf_pos == pos)[0][0]
                aims_idx = np.where(chr_aims_pos == pos)[0][0]
                vcf_indices.append(vcf_idx)
                aims_indices.append(aims_idx)

            # Extract data for overlapping AIMs
            selected_genotypes = chr_vcf_geno[vcf_indices]
            selected_ref = chr_vcf_ref[vcf_indices]
            selected_alt = chr_vcf_alt[vcf_indices]
            selected_anc1 = chr_aims_anc1[aims_indices]
            selected_anc2 = chr_aims_anc2[aims_indices]
            selected_positions = sorted(common_pos)

            # Store results for this chromosome
            ancestry_results[chr_name] = {
                'positions': selected_positions,
                'genotypes': selected_genotypes,
                'ref_alleles': selected_ref,
                'alt_alleles': selected_alt,
                'anc1_alleles': selected_anc1,  # gambiae or gambcolu
                'anc2_alleles': selected_anc2,  # coluzzii or arab
                'vcf_indices': vcf_indices,
                'aims_indices': aims_indices
            }

            print(f"  {len(common_pos)} AIMs processed")

        return ancestry_results

    # Run analysis for both datasets
    gambcolu_ancestry = analyze_ancestry_for_dataset(gambcolu_overlap_results, aims_gambcolu, "gambcolu")

    if arab_overlap_results and aims_arab is not None:
        arab_ancestry = analyze_ancestry_for_dataset(arab_overlap_results, aims_arab, "arabiensis")
    else:
        arab_ancestry = {}

    print(f"\nAncestry analysis complete!")
    print(f"Gambiae/Coluzzii: {sum([len(data['positions']) for data in gambcolu_ancestry.values()])} total AIMs")
    if arab_ancestry:
        print(f"Arabiensis: {sum([len(data['positions']) for data in arab_ancestry.values()])} total AIMs")
    return


@app.cell
def _(
    aims_arab,
    aims_gambcolu,
    arab_overlap_results,
    gambcolu_overlap_results,
    np,
    target_contigs,
    vcf_chromosomes,
    vcf_dataset,
    vcf_genotypes,
    vcf_positions,
    vcf_sample_names,
):
    def create_aims_genotype_matrix(overlap_results, aims_dataset, dataset_name):
        """Create genotype matrix for AIM visualization"""

        all_positions = []
        all_genotype_calls = []
        chromosome_boundaries = {}

        for chr_name in target_contigs:
            if chr_name not in overlap_results:
                continue

            print(f"Processing {dataset_name} chromosome {chr_name}...")

            # Get VCF data for this chromosome
            chr_mask = vcf_chromosomes == chr_name
            chr_vcf_pos = vcf_positions[chr_mask]
            chr_vcf_geno = vcf_genotypes[chr_mask]
            chr_vcf_ref = vcf_dataset['variants/REF'][chr_mask]
            chr_vcf_alt = vcf_dataset['variants/ALT'][chr_mask]

            # Get AIMs data
            chr_aims_pos = aims_dataset[chr_name]['POS'][:]
            if dataset_name == "gambcolu":
                chr_aims_anc1 = aims_dataset[chr_name]['gamb_allele'][:]  # gambiae
                chr_aims_anc2 = aims_dataset[chr_name]['colu_allele'][:]  # coluzzii
            else:  # arabiensis
                chr_aims_anc1 = aims_dataset[chr_name]['gambcolu_allele'][:]  # gambiae+coluzzii
                chr_aims_anc2 = aims_dataset[chr_name]['arab_allele'][:]     # arabiensis

            # Find matching positions
            common_pos = overlap_results[chr_name]['common_positions']
            sorted_common_pos = sorted(common_pos)

            # Store chromosome boundary
            start_idx = len(all_positions)
            chromosome_boundaries[chr_name] = start_idx

            # Get indices for matching positions
            chr_genotype_calls = []

            for pos in sorted_common_pos:
                vcf_idx = np.where(chr_vcf_pos == pos)[0][0]
                aims_idx = np.where(chr_aims_pos == pos)[0][0]

                # Get genotype data
                genotype = chr_vcf_geno[vcf_idx]
                ref_allele = chr_vcf_ref[vcf_idx]
                alt_alleles = chr_vcf_alt[vcf_idx]
                anc1_allele = chr_aims_anc1[aims_idx]
                anc2_allele = chr_aims_anc2[aims_idx]

                # Convert genotypes to ancestry calls for each sample
                pos_calls = []
                for sample_idx in range(len(vcf_sample_names)):
                    sample_gt = genotype[sample_idx]

                    # Convert genotype indices to nucleotides
                    alleles = []
                    for allele_idx in range(2):  # diploid
                        gt_val = sample_gt[allele_idx]
                        if gt_val == -1:  # missing
                            alleles.append('N')
                        elif gt_val == 0:  # reference
                            alleles.append(ref_allele)
                        else:  # alternative
                            alt_idx = gt_val - 1
                            if alt_idx < len(alt_alleles) and alt_alleles[alt_idx]:
                                alleles.append(alt_alleles[alt_idx])
                            else:
                                alleles.append('N')

                    # Classify ancestry based on alleles
                    anc1_count = sum(1 for a in alleles if a == anc1_allele)
                    anc2_count = sum(1 for a in alleles if a == anc2_allele)

                    if anc1_count == 2:
                        if dataset_name == "gambcolu":
                            call = 2  # gamb/gamb (blue)
                        else:
                            call = 2  # gambcolu/gambcolu (blue)
                    elif anc2_count == 2:
                        if dataset_name == "gambcolu":
                            call = 0  # colu/colu (red)
                        else:
                            call = 0  # arab/arab (red)
                    elif anc1_count == 1 and anc2_count == 1:
                        call = 1  # heterozygous (yellow)
                    else:
                        call = 3  # missing/other (white)

                    pos_calls.append(call)

                chr_genotype_calls.append(pos_calls)
                all_positions.append(pos)

            all_genotype_calls.extend(chr_genotype_calls)
            print(f"  {len(sorted_common_pos)} AIMs processed")

        # Convert to numpy array for plotting
        genotype_matrix = np.array(all_genotype_calls).T  # Transpose: samples x positions

        return genotype_matrix, all_positions, chromosome_boundaries

    # Create genotype matrices
    print("Creating gambiae/coluzzii genotype matrix...")
    gambcolu_matrix, gambcolu_positions, gambcolu_boundaries = create_aims_genotype_matrix(
        gambcolu_overlap_results, aims_gambcolu, "gambcolu")

    if arab_overlap_results and aims_arab is not None:
        print("Creating arabiensis genotype matrix...")
        arab_matrix, arab_positions, arab_boundaries = create_aims_genotype_matrix(
            arab_overlap_results, aims_arab, "arabiensis")
    else:
        arab_matrix = None

    print(f"\nMatrices created!")
    print(f"Gambiae/Coluzzii matrix shape: {gambcolu_matrix.shape} (samples x AIMs)")
    if arab_matrix is not None:
        print(f"Arabiensis matrix shape: {arab_matrix.shape} (samples x AIMs)")
    return (
        arab_boundaries,
        arab_matrix,
        arab_positions,
        gambcolu_boundaries,
        gambcolu_matrix,
        gambcolu_positions,
    )


@app.cell
def _(
    arab_boundaries,
    arab_matrix,
    arab_positions,
    gambcolu_boundaries,
    gambcolu_matrix,
    gambcolu_positions,
    plt,
    vcf_sample_names,
):
    import matplotlib.patches as patches

    def plot_aims_heatmap(genotype_matrix, positions, boundaries, dataset_name, sample_names):
        """Create AIM heatmap similar to your example image"""

        fig, ax = plt.subplots(figsize=(20, 8))

        # Define distinct colors for genotype calls
        if dataset_name == "gambcolu":
            colors = ['#D2691E', '#FFD700', '#4169E1', 'white']  # brown, yellow, royal blue, white
            labels = ['colu/colu', 'gamb/colu', 'gamb/gamb', 'missing']
        else:
            colors = ['#228B22', '#FFD700', '#9370DB', 'white']  # forest green, yellow, medium purple, white
            labels = ['arab/arab', 'gambcolu/arab', 'gambcolu/gambcolu', 'missing']

        # Create custom colormap
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)

        # Plot the heatmap
        im = ax.imshow(genotype_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=3)

        # Add chromosome boundaries
        boundary_positions = []
        chr_centers = []

        for i, (chr_name, start_idx) in enumerate(boundaries.items()):
            if i == 0:
                # First chromosome starts at 0
                boundary_positions.append(start_idx - 0.5)
            else:
                # Add vertical line at chromosome boundary
                ax.axvline(x=start_idx - 0.5, color='black', linewidth=3)
                boundary_positions.append(start_idx - 0.5)

            # Calculate center position for chromosome label
            if i < len(boundaries) - 1:
                next_chr = list(boundaries.keys())[i + 1]
                next_start = boundaries[next_chr]
                center = (start_idx + next_start) / 2
            else:
                center = (start_idx + len(positions)) / 2

            chr_centers.append((center, chr_name))

        # Add chromosome labels at the top
        for center, chr_name in chr_centers:
            ax.text(center, -0.8, chr_name, ha='center', va='bottom', 
                    fontsize=16, fontweight='bold')

        # Customize axes
        ax.set_xlabel('Variants', fontsize=14)
        ax.set_ylabel('Samples', fontsize=14)
        ax.set_title(f'AIM Genotypes - {dataset_name.title()} Analysis', fontsize=16, pad=20)

        # Set sample names on y-axis
        ax.set_yticks(range(len(sample_names)))
        ax.set_yticklabels(sample_names, fontsize=12)

        # Remove x-axis ticks (too many positions)
        ax.set_xticks([])

        # Add legend with larger font
        legend_elements = [patches.Patch(color=colors[i], label=labels[i]) for i in range(4)]
        ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), 
                  title='AIM genotype', fontsize=14, title_fontsize=14)

        # Adjust layout
        plt.tight_layout()

        return fig

    # Plot gambiae/coluzzii heatmap
    print("Creating gambiae/coluzzii AIM heatmap...")
    fig1 = plot_aims_heatmap(gambcolu_matrix, gambcolu_positions, gambcolu_boundaries, 
                            "gambcolu", vcf_sample_names)
    plt.show()



    # Plot arabiensis heatmap if available
    if arab_matrix is not None:
        print("\nCreating arabiensis AIM heatmap...")
        fig2 = plot_aims_heatmap(arab_matrix, arab_positions, arab_boundaries, 
                                "arabiensis", vcf_sample_names)
        plt.show()



    print("\nHeatmap visualization complete!")
    return fig1, fig2


@app.cell
def _(fig1, fig2):

    figures_info = [
        (fig1, "gambiae_coluzzii_AIM_heatmap.pdf"),
    ]

    if 'fig2' in locals() and fig2 is not None:
        figures_info.append((fig2, "arabiensis_AIM_heatmap.pdf"))

    # Save each figure
    for figure, filename in figures_info:
        figure.savefig(filename, format='pdf', dpi=300, bbox_inches='tight', 
                      facecolor='white', edgecolor='none')
        print(f"✅ Saved {filename}")

    print("🎉 All figures saved as PDF!")
    return


@app.cell
def _(mo):
    mo.md(r"""# Karyotyping plotting inversion markers""")
    return


@app.cell
def _(np, pd, vcf_chromosomes):
    # Load the inversion tags CSV file
    print("Loading inversion tags CSV file...")
    try:
        # Replace with your actual path
        inversion_tags_path = "/Users/talal.alyazeedi/Documents/LSTM_Postdoc/A.gambiae_BSA/rna-seq-pop/resources/karyotype_tag_snps.csv"  
    
        # Load the CSV
        df_tag_snps = pd.read_csv(inversion_tags_path)
    
        print("Inversion tags loaded successfully!")
        print(f"Shape: {df_tag_snps.shape}")
        print(f"Columns: {list(df_tag_snps.columns)}")
    
        # Display first few rows
        print("\nFirst 5 rows:")
        print(df_tag_snps.head())
    
        # Show available inversions
        if 'inversion' in df_tag_snps.columns:
            unique_inversions = df_tag_snps['inversion'].unique()
            print(f"\nAvailable inversions: {unique_inversions}")
        
            # Count tag SNPs per inversion
            print("\nTag SNPs per inversion:")
            inversion_counts = df_tag_snps['inversion'].value_counts()
            for inversion_name, count in inversion_counts.items():
                print(f"  {inversion_name}: {count} tag SNPs")
    
        # Check if chromosomes/contigs match between VCF and tag SNPs
        if 'contig' in df_tag_snps.columns or 'chromosome' in df_tag_snps.columns:
            contig_col = 'contig' if 'contig' in df_tag_snps.columns else 'chromosome'
            tag_contigs = df_tag_snps[contig_col].unique()
            print(f"\nContigs in tag SNPs: {tag_contigs}")
        
            # Check overlap with VCF contigs
            if 'vcf_chromosomes' in locals():
                vcf_contigs = np.unique(vcf_chromosomes)
                overlap = np.intersect1d(tag_contigs, vcf_contigs)
                print(f"Overlapping contigs: {overlap}")
    
        # Show position range
        if 'position' in df_tag_snps.columns:
            min_pos = df_tag_snps['position'].min()
            max_pos = df_tag_snps['position'].max()
            print(f"\nPosition range: {min_pos:,} - {max_pos:,}")
        
            # Show position range per contig
            if 'contig' in df_tag_snps.columns or 'chromosome' in df_tag_snps.columns:
                contig_col = 'contig' if 'contig' in df_tag_snps.columns else 'chromosome'
                print("\nPosition ranges per contig:")
                for chr_name in df_tag_snps[contig_col].unique():
                    contig_data = df_tag_snps[df_tag_snps[contig_col] == chr_name]
                    min_pos = contig_data['position'].min()
                    max_pos = contig_data['position'].max()
                    print(f"  {chr_name}: {min_pos:,} - {max_pos:,}")
    
    except FileNotFoundError:
        print(f"Error: File not found at {inversion_tags_path}")
        print("Please update the 'inversion_tags_path' variable with the correct path")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        df_tag_snps = None
    return (df_tag_snps,)


@app.cell
def _(df_tag_snps, np, pd, vcf_dataset):
    # Karyotype Analysis Function
    def analyze_karyotypes(vcf_data, tag_snps_df, selected_inversions=None, ploidy=2):
        """
        Analyze karyotypes using tag SNPs from VCF data
        """
        if selected_inversions is None:
            selected_inversions = ['2La', '2Rb']  # Default inversions
    
        print(f"Analyzing inversions: {selected_inversions}")
    
        # Extract VCF data
        vcf_pos = vcf_data['variants/POS']
        vcf_chrom = vcf_data['variants/CHROM'] 
        vcf_alt = vcf_data['variants/ALT']
        vcf_gt = vcf_data['calldata/GT']
        vcf_samples = vcf_data['samples']
    
        print(f"Working with {len(vcf_samples)} samples and {len(vcf_pos):,} variants")
    
        results_list = []
    
        # Process each selected inversion
        for inv in selected_inversions:
            print(f"\nProcessing inversion: {inv}")
        
            # Get tag SNPs for this inversion
            inv_tags = tag_snps_df[tag_snps_df['inversion'] == inv].copy()
        
            if len(inv_tags) == 0:
                print(f"  No tag SNPs found for {inv}")
                continue
            
            inv_contig = inv_tags['contig'].iloc[0]
            print(f"  Found {len(inv_tags)} tag SNPs on contig {inv_contig}")
        
            # Filter VCF to this contig
            contig_mask = vcf_chrom == inv_contig
            contig_pos = vcf_pos[contig_mask]
            contig_alt = vcf_alt[contig_mask]
            contig_gt = vcf_gt[contig_mask]
        
            print(f"  Contig {inv_contig} has {len(contig_pos):,} variants in VCF")
        
            # Find matching positions between tag SNPs and VCF
            tag_positions = inv_tags['position'].values
            matching_indices = []
            matching_tag_indices = []
        
            for i, tag_pos in enumerate(tag_positions):
                vcf_matches = np.where(contig_pos == tag_pos)[0]
                if len(vcf_matches) > 0:
                    matching_indices.append(vcf_matches[0])
                    matching_tag_indices.append(i)
        
            if len(matching_indices) == 0:
                print(f"  No matching positions found for {inv}")
                continue
            
            print(f"  Found {len(matching_indices)} matching tag SNP positions")
        
            # Extract matching data
            match_gt = contig_gt[matching_indices]
            match_alt = contig_alt[matching_indices]
            match_tag_alts = inv_tags['alt_allele'].iloc[matching_tag_indices].values
        
            # Calculate alternative allele counts for each sample
            sample_scores = []
            sample_totals = []
        
            for sample_idx in range(len(vcf_samples)):
                total_alt_count = 0
                total_sites = 0
            
                for site_idx in range(len(matching_indices)):
                    gt_calls = match_gt[site_idx, sample_idx]
                    site_alts = match_alt[site_idx]
                    target_alt = match_tag_alts[site_idx]
                
                    # Skip if missing data
                    if np.any(gt_calls == -1):
                        continue
                
                    # Find which alt allele index matches our target
                    target_alt_idx = None
                    if isinstance(site_alts, (list, np.ndarray)):
                        for alt_idx, alt in enumerate(site_alts):
                            if alt == target_alt:
                                target_alt_idx = alt_idx + 1  # +1 because 0 is reference
                                break
                
                    if target_alt_idx is not None:
                        # Count how many alleles match the target
                        alt_count = np.sum(gt_calls == target_alt_idx)
                        total_alt_count += alt_count
                        total_sites += 1
            
                sample_scores.append(total_alt_count)
                sample_totals.append(total_sites)
        
            # Calculate frequencies
            mean_scores = np.array(sample_scores) / np.maximum(np.array(sample_totals), 1)
            frequencies = mean_scores / ploidy
        
            # Create results for this inversion
            inv_results = pd.DataFrame({
                'sample_id': vcf_samples,
                'inversion': inv,
                f'karyotype_{inv}_mean': np.round(mean_scores, 2),
                f'karyotype_{inv}_freq': np.round(frequencies, 2),
                f'karyotype_{inv}_n_sites': sample_totals
            })
        
            print(f"  Results: freq range {frequencies.min():.2f} - {frequencies.max():.2f}")
            results_list.append(inv_results.set_index('sample_id'))
    
        if not results_list:
            print("No results generated!")
            return None
    
        # Combine all results
        karyotype_results = pd.concat(results_list, axis=1)
        print(f"\nFinal results shape: {karyotype_results.shape}")
    
        return karyotype_results

    # Run the analysis
    print("Starting karyotype analysis...")
    selected_inversions_to_analyze = ['2La', '2Rb']  # You can modify this list

    karyotype_df = analyze_karyotypes(
        vcf_dataset, 
        df_tag_snps, 
        selected_inversions=selected_inversions_to_analyze,
        ploidy=10
    )

    if karyotype_df is not None:
        print("\n" + "="*50)
        print("KARYOTYPE ANALYSIS RESULTS")
        print("="*50)
        print(f"Shape: {karyotype_df.shape}")
        print(f"Columns: {list(karyotype_df.columns)}")
        print("\nFirst 10 samples:")
        print(karyotype_df.head(10))
    
        # Summary statistics
        freq_cols = [col for col in karyotype_df.columns if 'freq' in col]
        if freq_cols:
            print(f"\nFrequency summary:")
            print(karyotype_df[freq_cols].describe())
    return (karyotype_df,)


@app.cell
def _(karyotype_df, plt, sns):
    def create_karyotype_figures_for_saving(karyo_df):
        """Create karyotype figures and return them for PDF saving"""
    
        # Extract frequency columns for visualization
        freq_columns = [col for col in karyo_df.columns if 'freq' in col]
        df_freq = karyo_df[freq_columns].copy()
    
        # Clean up column names for better display
        df_freq_display = df_freq.copy()
        df_freq_display.columns = [col.replace('karyotype_', '').replace('_freq', '') for col in df_freq_display.columns]
    
        figures = []
    
        # Figure 1: Main heatmap
        fig1 = plt.figure(figsize=(14, 6))
        ax1 = sns.heatmap(
            df_freq_display.T,  # Transpose so samples are on x-axis
            annot=True, 
            cmap='OrRd', 
            vmin=0, 
            vmax=1,
            fmt='.3f',
            cbar_kws={'label': 'Inversion Frequency'},
            linewidths=0.5,
            linecolor='white'
        )
        plt.title('Anopheles gambiae Karyotype Frequencies', fontsize=16, pad=20)
        plt.xlabel('Samples', fontsize=12)
        plt.ylabel('Inversions', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.show()  # Show the figure
        figures.append((fig1, "karyotype_heatmap.pdf"))
    
        # Figure 2: Distribution plots
        fig2 = plt.figure(figsize=(12, 5))
    
        # Plot 1: Bar plot of frequencies by sample
        plt.subplot(1, 2, 1)
        x_positions = range(len(df_freq_display))
        width = 0.35
        for i, col in enumerate(df_freq_display.columns):
            plt.bar([x + i*width for x in x_positions], df_freq_display[col], width, 
                   label=col, alpha=0.8)
        plt.xlabel('Samples')
        plt.ylabel('Inversion Frequency')
        plt.title('Karyotype Frequencies by Sample')
        plt.xticks([x + width/2 for x in x_positions], df_freq_display.index, rotation=45, ha='right')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
        # Plot 2: Distribution of frequencies
        plt.subplot(1, 2, 2)
        for col in df_freq_display.columns:
            plt.hist(df_freq_display[col], bins=10, alpha=0.7, label=col, edgecolor='black')
        plt.xlabel('Inversion Frequency')
        plt.ylabel('Number of Samples')
        plt.title('Distribution of Inversion Frequencies')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()  # Show the figure
        figures.append((fig2, "karyotype_distribution.pdf"))
    
        return figures

    # Create and show figures, then assign to variables
    karyotype_figures = create_karyotype_figures_for_saving(karyotype_df)

    # Assign individual figures to variables for easy access
    karyotype_heatmap_fig = karyotype_figures[0][0]
    karyotype_distribution_fig = karyotype_figures[1][0]

    print(f"\nCreated and displayed {len(karyotype_figures)} figures")
    print("Variables created:")
    print("  • karyotype_heatmap_fig - Main heatmap figure")
    print("  • karyotype_distribution_fig - Distribution plots figure") 
    print("  • karyotype_figures - List of (figure, filename) tuples")
    print("\nReady for saving in next cell!")
    return karyotype_distribution_fig, karyotype_heatmap_fig


@app.cell
def _(karyotype_distribution_fig, karyotype_heatmap_fig, plt):
    figures_karyo = [
        (karyotype_heatmap_fig, "karyotype_heatmap.pdf"),
        (karyotype_distribution_fig, "karyotype_distribution.pdf")
    ]

    # Save each figure
    for figs, file_name in figures_karyo:
        figs.savefig(file_name, format='pdf', dpi=300, bbox_inches='tight', 
                      facecolor='white', edgecolor='none')
        print(f"✅ Saved {file_name}")

    print("🎉 All karyotype figures saved as PDF!")

    # Close figures to free memory
    for figs, file_name in figures_karyo:
        plt.close(figs)
    return


if __name__ == "__main__":
    app.run()
