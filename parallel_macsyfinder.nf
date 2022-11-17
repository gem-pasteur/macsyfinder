#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*************************
 * Default options
 *************************/

params['models'] = false
params['sequence-db'] = false
params['db-type'] = false
params['replicon-topoloy']= false
params['topology-file']= false
params['inter-gene-max-space']= false
params['min-mandatory-genes-required']= false
params['min-genes-required']= false
params['max-nb-genes']= false
params['multi-loci']= false
params['hmmer']= false
params['e-value-search']= false
params['no-cut-ga']= false
params['i-evalue-sel']= false
params['coverage-profile']= false
params['mandatory-weight']= false
params['accessory-weight']= false
params['exchangeable-weight']= false
params['redundancy-penalty']= false
params['out-of-cluster']= false
params['models-dir']= false
params['index-dir']= false
params['res-search-suffix']= false
params['res-extract-suffix']= false
params['profile-suffix']= false
params['worker']= false
params['outdir']= false
params['cfg-file']= false
params['debug'] = false

/****************************************
*  parameter not really used in the wf  *
*  used to improve UI                   *
*****************************************/
params.help = false

// the nextflow profile
params.profile = false

/****************************************
*             real parameters           *
*****************************************/

models = params.models ? " --models ${params['models']}" : ''
sequence_db = params['sequence-db'] ? "${params['sequence-db']}": ''
db_type = params['db-type']
replicon_topology = params['replicon-topoloy'] ? " --replicon-topology ${params['replicon-topoloy']}" : ''
topology_file = params['topology-file'] ? " --topology-file ${params['topology-file']}" : ''
inter_gene_max_space = params['inter-gene-max-space'] ? " --inter-gene-max-space ${params['inter-gene-max-space']}" : ''
min_mandatory_genes_required = params['min-mandatory-genes-required'] ? " --min-mandatory-genes-required ${params['min-mandatory-genes-required']}": ''
min_genes_required = params['min-genes-required'] ? " --min-genes-required ${params['min-genes-required']}" : ''
max_nb_genes = params['max-nb-genes'] ? " --max-nb-genes ${params['max-nb-genes']}" : ''
multi_loci = params['multi-loci'] ? " --multi-loci ${params['multi-loci']}" : ''
hmmer = params['hmmer'] ? " --hmmer ${params['hmmer']}" : ''
e_value_search = params['e-value-search'] ? " --e-value-search ${params['e-value-search']}" : ''
no_cut_ga = params['no-cut-ga'] ? ' --no-cut-ga' : ''
i_value_sel = params['i-evalue-sel'] ? " --i-evalue-sel ${params['i-evalue-sel']}" : ''
coverage_profile = params['coverage-profile'] ? " --coverage-profile ${params['coverage-profile']}" : ''
mandatory_weight = params['mandatory-weight'] ? " --mandatory-weight ${params['mandatory-weight']}" : ''
accessory_weight = params['accessory-weight'] ? " --accessory-weight ${params['accessory-weight']}" : ''
exchangeable_weight = params['exchangeable-weight'] ? " --exchangeable-weight ${params['exchangeable-weight']}" : ''
redundancy_penalty = params['redundancy-penalty'] ? " --redundancy-penalty ${params['redundancy-penalty']}" : ''
out_of_cluster = params['out-of-cluster'] ? " --out-of-cluster ${params['out-of-cluster']}" : ''
models_dir = params['models-dir'] ? " --models-dir ${file(params['models-dir'])}" : ''
index_dir = params['index-dir'] ? " --index-dir ${params['index-dir']}" : ''
res_search_suffix = params['res-search-suffix'] ? " --res-search-suffix ${params['res-search-suffix']}" : ''
res_extract_suffix = params['res-extract-suffix'] ? " --res-extract-suffix ${params['res-extract-suffix']}" : ''
profile_suffix = params['profile-suffix'] ? " --profile-suffix ${params['profile-suffix']}" : ''
worker = params['worker'] ? " --worker ${params['worker']}" : ''
cfg_file= params['cfg-file'] ? " --cfg-file ${params['cfg-file']}" : ''
debug = params.debug ? ' -vv' : ''


if (params.help){
   msg = '''
parallel_macsyfinder available options:

 --models
 --sequence-db
 --db-type
 --replicon-topology
 --topology-file
 --inter-gene-max-space
 --min-mandatory-genes-required
 --min-genes-required
 --max-nb-genes
 --multi-loci
 --hmmer
 --e-value-search
 --no-cut-ga
 --i-evalue-sel
 --coverage-profile
 --mandatory-weight
 --accessory-weight
 --exchangeable-weight
 --redundancy-penalty
 --out-of-cluster-weight
 --models-dir
 --index-dir
 --res-search-suffix
 --res-extract-suffix
 --profile-suffix
 --cfg-file
 --worker
 --outdir
Please refer to the MacSyFinder documentation (https://macsyfinder.readthedocs.io) for the meaning of each options.
'''
    println(msg)
    System.exit(0)
}

if (params.profile){
    throw new Exception("The macsyfinder option '--profile' does not exists. May be you want to use the nextflow option '-profile'.")
}

if (!params['db-type']){
    throw new Exception("The option '--db-type' is mandatory.")
} else if (!params['db-type'] in ['gembase', 'ordered_replicon', 'unordered']) {
    throw new Exception("The option '--db-type' accept only 'gembase', 'ordered_replicon', 'unordered' as value: get '${params['db-type']}'")
}

if (! params['sequence-db'] and ! params['cfg-file']){
    throw new Exception("The option '--sequence-db' is mandatory.");
}
if (! params['models'] and ! params['cfg-file']){
    throw new Exception("The option '--models' is mandatory.");
}

if ( params['db-type'] != 'gembase' && ! params['outdir'] ){
    throw new Exception("The option 'outdir' is mandatory if db-type = '${params['db-type']}''")
}

sequence_db = file(sequence_db)

if( params.outdir ){
    outdir = params.outdir
} else if ( params['db-type'] == 'gembase' ){
    outdir = "merged_macsyfinder_results_${sequence_db.baseName}"
} else {
    throw new Exception("The option 'outdir' is manadatory if db-type != 'gembase'")
}

/****************************************
 *           The workflow               *
 ****************************************/

// split the gembase file in replicon files
process split{

    input:
        path gembase
    output:
        path "split_gembase/*.fasta"
    script:
        """
        macsysplit --mute --outdir split_gembase ${gembase}
        """
}

// execute macsyfinder on each replicon file
process macsyfinder{

    input:
        path one_replicon
        val models
        val db_type
        val replicon_topology
        val topology_file
        val inter_gene_max_space
        val min_mandatory_genes_required
        val min_genes_required
        val max_nb_genes
        val multi_loci
        val hmmer
        val e_value_search
        val no_cut_ga
        val i_value_sel
        val coverage_profile
        val mandatory_weight
        val accessory_weight
        val exchangeable_weight
        val redundancy_penalty
        val out_of_cluster
        val models_dir
        val index_dir
        val res_search_suffix
        val res_extract_suffix
        val profile_suffix
        val cfg_file
        val debug
    output:
        path("macsyfinder-${one_replicon.baseName}")

    script:
        """
        macsyfinder --sequence-db ${one_replicon} --db-type ${db_type} ${replicon_topology}${topology_file}${models_dir}${models}${inter_gene_max_space}${min_mandatory_genes_required}${min_genes_required}${max_nb_genes}${multi_loci} \
${hmmer}${e_value_search} ${no_cut_ga}${i_value_sel} ${coverage_profile}${mandatory_weight}${accessory_weight}${exchangeable_weight}${redundancy_penalty}${out_of_cluster} \
${index_dir}${res_search_suffix}${res_extract_suffix}${profile_suffix} --worker ${task.cpus} --mute${debug} --out-dir macsyfinder-${one_replicon.baseName}
"""
}


// merge the results of macsyfinder on each replicon.
process merge_results{
    publishDir outdir, mode: 'copy'

    input:
        path all_results_dirs

    output:
        path "merged_*.tsv"
        path "merged_*.txt"
        path "macsy_merge_results.out"

    script:
        """
        macsymerge --mute ${all_results_dirs}
        """
}


workflow {

    if( params['db-type'] == 'gembase' ){
        gembase = Channel.fromPath(sequence_db)
        replicons = split(gembase)
    } else {
        if( params['sequence-db'].contains(',') ){
            paths = params['sequence-db'].tokenize(',')
        } else {
            paths = params['sequence-db']
        }
        replicons = Channel.fromPath(paths)
    }

    results_per_replicon = macsyfinder(replicons.flatten(), models, db_type, replicon_topology, topology_file,
    inter_gene_max_space, min_mandatory_genes_required, min_genes_required, max_nb_genes, multi_loci,
    hmmer, e_value_search, no_cut_ga, i_value_sel, coverage_profile,
    mandatory_weight, accessory_weight, exchangeable_weight, redundancy_penalty, out_of_cluster,
    models_dir, index_dir, res_search_suffix, res_extract_suffix,profile_suffix,
    cfg_file, debug)

    results = merge_results(results_per_replicon.toList())
}

workflow.onComplete {
    if ( workflow.success )
        println("\nDone!")
}

workflow.onError {
    println "Oops .. something went wrong"
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
