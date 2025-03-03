convert_to_tensor <- function(lst) {
  torch <- reticulate::import("torch")
  convert_df <- function(df) {
    df[] <- lapply(df, function(x) {
      if (is.factor(x)) as.numeric(as.character(x)) else x
    })
    # Ensure all data is numeric, replacing non-numeric with NA
    df[] <- lapply(df, function(x) {
      as.numeric(x)
    })
    torch$from_numpy(as.matrix(df))
  }
  
  convert_list <- function(x) {
    if (is.data.frame(x)) {
      return(convert_df(x))
    } else if (is.list(x)) {
      return(lapply(x, convert_list))
    } else {
      return(x)
    }
  }
  convert_list(lst)
}
