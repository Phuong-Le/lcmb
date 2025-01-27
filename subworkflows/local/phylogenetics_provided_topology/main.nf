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
    // assign mutation to provided topology
    mutToTree(
        phylogenetics_provided_topology_input_ch,
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
