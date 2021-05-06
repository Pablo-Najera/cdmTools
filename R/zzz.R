.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

StartWelcomeMessage <- function(){
  paste(c("==========================================================\n",
          "cdmTools Package",
          " [Version ", utils::packageDescription("cdmTools")$Version,
          "; ",utils::packageDescription("cdmTools")$Date, "]\n",
          "More information: https://github.com/Pablo-Najera/cdmTools\n",
          "==========================================================\n"),
        sep="")
}
