# This simple script tests R dependency installation for VDJtools
args    <- commandArgs(TRUE)

# iterate over dependencies
passed <- T
for (i in 1:length(args)) {
   if (!args[i] %in% rownames(installed.packages())) {
       print(paste("FAILED to install ", args[i]))
       passed <- F
   }
}

if (passed) {
   print("PASSED")
}

