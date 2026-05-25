// Helper functions for pipeline tests

class UTILS {

    public static def getAssertions = { Map args ->
        // Mandatory, as we always need an outdir
        def outdir = args.outdir

        // Get scenario and extract all properties dynamically
        def scenario = args.scenario ?: [:]

        // Pass down workflow for std capture
        def workflow = args.workflow

        // These strings are not stable and should be ignored
        def snapshot_ignore_list = [
            "Creating env using",
            "Downloading plugin",
            "Got an interrupted  exception while taking agent result",
            "Pulling Singularity image",
            "Staging foreign file",
            "Unable to resume cached task",
            "Unable to stage foreign file",
        ]

        // stable_name: All files + folders in ${outdir}/ with a stable name
        def stable_name = getAllFilesFromPath(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_content: All files in ${outdir}/ with stable content
        def stable_content = getAllFilesFromPath(outdir, ignoreFile: 'tests/.nftignore', ignore: [scenario.ignoreFiles])

        // bam_files: All bam files
        def bam_files = getAllFilesFromPath(outdir, include: ['**/*.bam'], ignore: [scenario.ignoreFiles])
        // recal_bam_files: All recalibrated bam files
        def recal_bam_files = getAllFilesFromPath(outdir, include: ['**/*.recal.bam']) - bam_files
        // cram_files: All cram files
        def cram_files = getAllFilesFromPath(outdir, include: ['**/*.cram'], ignore: [scenario.ignoreFiles])
        // Fasta file for cram verification with nft-bam
        def fasta_base = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
        def fasta = fasta_base + 'genomics/homo_sapiens/genome/genome.fasta'
        // vcf_files: All vcf files
        def vcf_files = getAllFilesFromPath(outdir, include: ['**/*.vcf{,.gz}'], ignore: [scenario.ignoreFiles])

        def assertion = []
        // getAllFilesFromPath returns relative paths (strings), so this resolves to an absolute path
        def absolutePath = { file -> file.toString().startsWith('/') ? file.toString() : "${outdir}/${file}" }

        if (!scenario.failure) {
            assertion.add(workflow.trace.succeeded().size())
            assertion.add(removeFromYamlMap("${outdir}/pipeline_info/nf_core_rnavar_software_mqc_versions.yml", "Workflow")?: 'No versions')
        }

        // At least always pipeline_info/ is created and stable
        assertion.add(stable_name)

        if (!scenario.stub) {
            assertion.add(stable_content.isEmpty() ? 'No stable content' : stable_content.collect { file -> path(absolutePath(file)) })
            assertion.add(bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.tokenize('/').last() + ":md5," + bam(absolutePath(file)).readsMD5 })
            assertion.add(recal_bam_files.isEmpty() ? 'No unstable recal BAM files' : recal_bam_files.collect { file -> file.tokenize('/').last() + ":stats" + bam(absolutePath(file)).getStatistics() })
            assertion.add(cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.tokenize('/').last() + ":md5," + cram(absolutePath(file), fasta).readsMD5 })
            assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.tokenize('/').last() + ":md5," + path(absolutePath(file)).vcf.variantsMD5 })
        }

        // If we have a snapshot options in scenario then we allow to capture either stderr, stdout or both
        // With options to include specific stings
        def workflow_std = []
        // Otherwise, we always capture stdout and stderr for any WARN message
        // Both have additional possibilities to ignore some strings
        def filter_args = [ignore: snapshot_ignore_list + (scenario.snapshot_ignore ?: [])]

        workflow_std = workflow.stderr + workflow.stdout
        filter_args.include = ["WARN"]

        assertion.add(filterNextflowOutput(workflow_std, filter_args) ?: "No warnings")

        if (scenario.snapshot) {
            workflow_std = scenario.snapshot.split(',')
                .findAll { it in ['stderr', 'stdout'] }
                .collect { workflow."$it" }
                .flatten()

            filter_args.remove('include')

            if (scenario.snapshot_include) {
                filter_args.include = [scenario.snapshot_include]
            }

            assertion.add(filterNextflowOutput(workflow_std, filter_args) ?: "No content")
        }

        return assertion
    }

    public static def getTest = { scenario ->
        // This function returns a closure that will be used to run the test and the assertion
        // It will create tags or options based on the scenario

        return {
            // If the test is for a gpu, we add the gpu tag
            // Otherwise, we add the cpu tag
            // If the tests has no conda incompatibilities
            // then we append "_conda" to the cpu/gpu tag
            // If the test is for a stub, we add options -stub
            // And we append "_stub" to the cpu/gpu tag

            // All options should be:
            // gpu (this is the default for gpu)
            // cpu (this is the default for tests without conda)
            // gpu_conda (this should never happen)
            // cpu_conda (this is the default for tests with conda compatibility)
            // gpu_stub
            // cpu_stub
            // gpu_conda_stub (this should never happen)
            // cpu_conda_stub

            tag "pipeline"
            tag "pipeline_rnavar"

            options "-output-dir ${outputDir}${scenario.stub ? ' -stub' : ''}"

            if (scenario.gpu) {
                tag "gpu${!scenario.no_conda ? '_conda' : ''}${scenario.stub ? '_stub' : ''}"
            }

            if (!scenario.gpu) {
                tag "cpu${!scenario.no_conda ? '_conda' : ''}${scenario.stub ? '_stub' : ''}"
            }

            // If a tag is provided, add it to the test
            if (scenario.tag) {
                tag scenario.tag
            }

            when {
                params {
                    // Mandatory, as we always need an outdir
                    outdir = "${outputDir}"
                    // Apply scenario-specific params
                    scenario.params.each { key, value ->
                        delegate."$key" = value
                    }
                }
            }

            then {
                // Assert failure/success, and fails early so we don't pollute console with massive diffs
                if (scenario.failure) {
                    assert workflow.failed
                } else {
                    assert workflow.success
                }
                assertAll(
                    { assert snapshot(
                        // All assertions based on the scenario
                        *UTILS.getAssertions(
                            outdir: params.outdir,
                            scenario: scenario,
                            workflow: workflow
                        )
                    ).match() }
                )
            }
            cleanup {
                if (System.getenv('NFT_CLEANUP')) {
                    println ""
                    println "CLEANUP"
                    println "Set NFT_CLEANUP to false to disable."
                    println "The following folder will be deleted:"
                    println "- ${launchDir}"

                    new File("${launchDir}").deleteDir()
                }
            }
        }
    }
}
