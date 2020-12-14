#' SHOW model class.
#'
#' @export
show_model <- R6::R6Class(
  classname = "show_model",
  inherit = psd_model,

  private = list(
    Temp_ = NULL, # temperature
    n_phi_ = 3,
    tmb_DLL_ = "realPSD_TMBExports",
    tmb_model_ = "SHOW_log"
  ),

  active = list(
    #' @field Temp Temperature of the system.
    Temp = function() private$Temp_
  )

)
