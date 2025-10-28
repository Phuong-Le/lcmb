include { rmPolyclonal } from "$projectDir/modules/local/rm_polyclonal"
include { mutToTree } from "$projectDir/modules/local/mut_to_tree"
include { matrixGeneratorBranches } from "$projectDir/modules/local/matrix_generator_branches"
include { concatMatrices } from "$projectDir/modules/local/concat_matrices"
include { spectraPlottingBranches } from "$projectDir/modules/local/spectra_plotting_branches"


workflow PHYLOGENETICS_PROVIDED_TREE_TOPOLOGY {
    take:
    phylogenetics_provided_topology_input_ch
    outdir_basename
    sigprofiler_genome

    main:
    topology = phylogenetics_provided_topology_input_ch
    .map { pdid, nr_path, nv_path, genotype_bin_path, clonality, topology_file ->
        tuple(pdid, topology_file)
    }

    phylogenetic_input_ch = phylogenetics_provided_topology_input_ch
    .map { pdid, nr_path, nv_path, genotype_bin_path, clonality, topology_file ->
        tuple(pdid, nr_path, nv_path, genotype_bin_path, clonality)
    }

    // remove polyclonal samples
    rmPolyclonal(
        phylogenetic_input_ch,
        outdir_basename
    )

    // assign mutation to provided topology
    // only keep trees with >2 samples
    mutToTree(
        rmPolyclonal.out.phylogenetic_input
        .combine( topology, by: 0 )
        .filter { it -> it[3].readLines().first().split(' ').size() > 2 },
        outdir_basename
    )

    // generate mutation matrix for the branches by SigProfilerMatrixGenerator
    matrixGeneratorBranches(
        mutToTree.out.muts_on_tree,
        outdir_basename,
        sigprofiler_genome
        )

    concatMatrices(
        matrixGeneratorBranches.out.toList(),
        outdir_basename
        )

    // // plotting
    // spectraPlottingBranches(
    //     concatMatrices.out,
    //     outdir_basename
    //     )
}
