#' Document data of this package

#' Circ
#' 
#' CircRNA expression count table. This is come from the running result of DCC.
#' @name Circ
#' @docType data
#' @format A data frame with 153 rows and 21 variables.
#' \describe{
#'  \item{Chr}{Chromosomes of CircRNAs}
#'  \item{Start}{Start positions of CircRNAs}
#'  \item{End}{End positions of CircRNAs}
#'  \item{Vrg_F_1_1}{Drosophila virgin female, eclosion + 1 day, heads, replicate 1}
#'  \item{Vrg_F_1_2}{Drosophila virgin female, eclosion + 1 day, heads, replicate 2}
#'  \item{Mat_F_1_1}{Drosophila mated female, eclosion + 1 day, heads, replicate 1}
#'  \item{Mat_F_1_2}{Drosophila mated female, eclosion + 1 day, heads, replicate 2}
#'  \item{Mat_M_1_1}{Drosophila mated male, eclosion + 1 day, heads, replicate 1}
#'  \item{Mat_M_1_2}{Drosophila mated male, eclosion + 1 day, heads, replicate 2}
#'  \item{Vrg_F_4_1}{Drosophila virgin female, eclosion + 4 day, heads, replicate 1}
#'  \item{Vrg_F_4_2}{Drosophila virgin female, eclosion + 4 day, heads, replicate 2}
#'  \item{Mat_F_4_1}{Drosophila mated female, eclosion + 4 day, heads, replicate 1}
#'  \item{Mat_F_4_2}{Drosophila mated female, eclosion + 4 day, heads, replicate 2}
#'  \item{Mat_M_4_1}{Drosophila mated male, eclosion + 4 day, heads, replicate 1}
#'  \item{Mat_M_4_2}{Drosophila mated male, eclosion + 4 day, heads, replicate 2}
#'  \item{Vrg_F_20_1}{Drosophila virgin female, eclosion + 20 day, heads, replicate 1}
#'  \item{Vrg_F_20_2}{Drosophila virgin female, eclosion + 20 day, heads, replicate 2}
#'  \item{Mat_F_20_1}{Drosophila mated female, eclosion + 20 day, heads, replicate 1}
#'  \item{Mat_F_20_2}{Drosophila mated female, eclosion + 20 day, heads, replicate 2}
#'  \item{Mat_M_20_1}{Drosophila mated male, eclosion + 20 day, heads, replicate 1}
#'  \item{Mat_M_20_2}{Drosophila mated male, eclosion + 20 day, heads, replicate 2}
#' }
#' @source Original data from Westholm et al., 2014
#' 
NULL


#' Linear
#' 
#' CircRNA host gene expression count table. This is come from the running result of DCC.
#' @name Linear
#' @docType data
#' @format A data frame with 153 rows and 21 variables.
#' \describe{
#'  \item{Chr}{Chromosomes of CircRNAs}
#'  \item{Start}{Start positions of CircRNAs}
#'  \item{End}{End positions of CircRNAs}
#'  \item{Vrg_F_1_1}{Drosophila virgin female, eclosion + 1 day, heads, replicate 1}
#'  \item{Vrg_F_1_2}{Drosophila virgin female, eclosion + 1 day, heads, replicate 2}
#'  \item{Mat_F_1_1}{Drosophila mated female, eclosion + 1 day, heads, replicate 1}
#'  \item{Mat_F_1_2}{Drosophila mated female, eclosion + 1 day, heads, replicate 2}
#'  \item{Mat_M_1_1}{Drosophila mated male, eclosion + 1 day, heads, replicate 1}
#'  \item{Mat_M_1_2}{Drosophila mated male, eclosion + 1 day, heads, replicate 2}
#'  \item{Vrg_F_4_1}{Drosophila virgin female, eclosion + 4 day, heads, replicate 1}
#'  \item{Vrg_F_4_2}{Drosophila virgin female, eclosion + 4 day, heads, replicate 2}
#'  \item{Mat_F_4_1}{Drosophila mated female, eclosion + 4 day, heads, replicate 1}
#'  \item{Mat_F_4_2}{Drosophila mated female, eclosion + 4 day, heads, replicate 2}
#'  \item{Mat_M_4_1}{Drosophila mated male, eclosion + 4 day, heads, replicate 1}
#'  \item{Mat_M_4_2}{Drosophila mated male, eclosion + 4 day, heads, replicate 2}
#'  \item{Vrg_F_20_1}{Drosophila virgin female, eclosion + 20 day, heads, replicate 1}
#'  \item{Vrg_F_20_2}{Drosophila virgin female, eclosion + 20 day, heads, replicate 2}
#'  \item{Mat_F_20_1}{Drosophila mated female, eclosion + 20 day, heads, replicate 1}
#'  \item{Mat_F_20_2}{Drosophila mated female, eclosion + 20 day, heads, replicate 2}
#'  \item{Mat_M_20_1}{Drosophila mated male, eclosion + 20 day, heads, replicate 1}
#'  \item{Mat_M_20_2}{Drosophila mated male, eclosion + 20 day, heads, replicate 2}
#' }
#' @source Original data from Westholm et al., 2014
#' 
NULL
#' Coordinates
#' @name Coordinates
#' @docType data
#' @format A data frame with 153 rows and 6 variables.
#' \describe{
#'  \item{Chr}{Chromosomes of CircRNAs}
#'  \item{Start}{Start positions of CircRNAs}
#'  \item{End}{End positions of CircRNAs}
#'  \item{gene}{gene name}
#'  \item{junctiontype}{junction type. 1 and 2 are GT-AG and reverse complementary. 0 is non canonical junction}
#'  \item{strand}{strand of circRNA}
#' }
NULL