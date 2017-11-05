workflow stats {
  call run
  
  output {
    File output_stats = run.output_stats
  }
}

task run {
  String assembly
  String outdir
  String stat_sh
  String asm_name = basename(assembly)
  command {
    ${stat_sh} ${assembly} > ${asm_name}.stats
  }
  output {
    File output_stats = "${asm_name}.stats"
  }
  runtime {
    memory: "1 GB"
    cpu: "1"
    jobname: "run_stats"
  }
}
