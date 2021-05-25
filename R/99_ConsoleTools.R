
.msg          <- function(...) cat(crayon::bold(crayon::white(paste0(...))))
.msg_alt      <- function(...) cat(crayon::italic(crayon::bold(paste0(...))))
.msg_val      <- function(...) cat(crayon::italic(crayon::bgBlack(crayon::white(paste0(...)))))
.msg_alt_good <- function(...) cat(crayon::green(crayon::italic(crayon::bold(paste0(...)))))
.msg_alt_bad  <- function(...) cat(crayon::red(crayon::italic(crayon::bold(paste0(...)))))
.msg_name     <- function(...) cat(crayon::bgWhite(crayon::bold(crayon::black(paste0(...)))))
.msg_lite     <- function(...) cat(crayon::italic(crayon::bgBlack(crayon::white(paste0(...)))))
.msg_hpc      <- function(...) cat(crayon::bgBlue(crayon::white$bold(paste0(...))))
.msg_python   <- function(...) cat(crayon::bgYellow(crayon::black$bold(paste0(...))))
  