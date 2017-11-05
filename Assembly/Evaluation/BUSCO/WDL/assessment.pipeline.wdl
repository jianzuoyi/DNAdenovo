import "busco.wdl"
import "stats.wdl"

workflow Assessment_assembly {
  String assembly  
  
  call stats.run as stats { input: assembly = assembly }
  call busco.run as busco { input: assembly = assembly }
}

