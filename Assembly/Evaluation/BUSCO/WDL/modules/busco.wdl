workflow busco {
  call run
}

task run {
  String assembly
  String busco_config_file
  String augustus_bin
  String augustus_script
  String augustus_config_path
  String OrthoDB
  String runBUSCO
  String outbase
  String python
  String cpu
  String memory_gb

  command {
    export BUSCO_CONFIG_FILE=${busco_config_file}
    export PATH=${augustus_bin}:${augustus_script}:$PATH
    export AUGUSTUS_CONFIG_PATH=${augustus_config_path}
    python ${runBUSCO} -i ${assembly} -o ${outbase} -l ${OrthoDB} -m genome -c ${cpu} -f -r
  }
  runtime {
    memory: "${memory_gb} GB"
    cpu: "${cpu}"
    jobname: "run_busco"
  }
}
