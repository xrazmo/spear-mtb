nextflow.enable.dsl=2

include {MAP_READS as mpr} from "$baseDir/workflows/clockwork/main"
include {REMOVE_CONTAM as rmc} from "$baseDir/workflows/clockwork/main"
include {VARIANT_CALL as vrc} from "$baseDir/workflows/clockwork/main"
include {PREDICT_DST as prd} from "$baseDir/workflows/clockwork/main"

include {TBPROFILER_PROFILE as tbp} from "$baseDir/workflows/tbprofiler/profile/main"

params.ref_dir = "${baseDir}/ref_db"
params.input_gz = ""
params.prefix = ""
params.input_catalog = "${baseDir}/catalogues"

def vcf_dir = "${baseDir}/out_vcf"
def out_dir = "${baseDir}/out_csv"
def h37Rv_dir = "${params.ref_dir}/Ref.H37Rv";

workflow{

   reads_ch = Channel.fromFilePairs("${params.input_gz}/*/*_{1,2}.fastq.gz").map{it->[[id:it[0]],it[1]]};
   ref_fa = Channel.fromPath("${params.ref_dir}/Ref.remove_contam/*.fa");

   // Running TBPROFILER
   tbp(reads_ch)

   // Running CRyPTIC workflow
   ref_tsv = Channel.fromPath("${params.ref_dir}/Ref.remove_contam/*.tsv");

   mpr(reads_ch.combine(ref_fa)) 
   rmc(mpr.out.sam.combine(ref_tsv))
   vrc(rmc.out.reads,h37Rv_dir)
   
   file(vcf_dir).mkdir()
   vrc.out.final_vcf.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)

    catalog_ch = Channel.fromPath("${params.input_catalog}/*/*.csv")
    refpkl_ch = Channel.fromPath("${params.input_catalog}/*/*.gz")

    ch = vrc.out.final_vcf.combine(catalog_ch).combine(refpkl_ch)
   
    prd(ch)
    
    file(vcf_dir).mkdir()
    prd.out.effects.map{it->it[1]}.flatten()
        .collectFile(name:"merged.effects.${params.prefix}.csv",keepHeader:true,storeDir:out_dir)

    // Merging results    

} 