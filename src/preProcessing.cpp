#include <Rcpp.h>
#include <unordered_map>
#include "constants.h"
#include "utils.h"
using namespace Rcpp;

//' Rename individuals in a pedigree from 1 to N
//'
//' Functions in \bold{purgeR} require individuals to be named from 1 to N.
//' This takes a dataframe containing a pedigree, and rename individuals having
//' names in any format to that required by other functions in \pkg{purgeR}.
//' 
//' @template ped-arg
//' @param id A string naming the column with individual identities. It will be renamed to its default value 'id'.
//' @param dam A string naming the column with maternal identities. It will be renamed to its default value 'dam'.
//' @param sire A string naming the column with paternal identities. It will be renamed to its default value 'sire'.
//' @param keep_names A boolean value indicating whether the original identity values should be kept on a separate column (named 'names'), or not.
//' @return A dataframe with the pedigree's identities renamed.
// [[Rcpp::export]]
DataFrame rename(DataFrame ped,
                 std::string id = "id",
                 std::string dam = "dam",
                 std::string sire = "sire",
                 bool keep_names = false) {

  // Set input columns
	std::string id_name = id_col;
	std::string dam_name = dam_col;
	std::string sire_name = sire_col;
	if (id != id_col) id_name = id;
	if (dam != dam_col) dam_name = dam;
	if (sire != sire_col) sire_name = sire;

	// Rename
	Rcpp::CharacterVector str_id = ped[id_name];
	Rcpp::CharacterVector str_dam = ped[dam_name];
	Rcpp::CharacterVector str_sire = ped[sire_name];
	std::map<std::string, int> map_id;
	Rcpp::IntegerVector int_id;
	Rcpp::IntegerVector int_dam;
	Rcpp::IntegerVector int_sire;
	R_xlen_t nRows = ped.rows();
	for (int i(0); i<nRows; ++i) {
	  int id_i = i+1;
	  map_id[Rcpp::as<std::string>(str_id(i))] = id_i;
	  int_id.push_back(id_i);
	  int dam_i, sire_i;
	  dam_i = map_kv(map_id, Rcpp::as<std::string>(str_dam(i)), id_unk_int);
	  sire_i = map_kv(map_id, Rcpp::as<std::string>(str_sire(i)), id_unk_int);
	  int_dam.push_back(dam_i);
	  int_sire.push_back(sire_i);
	}

	Rcpp::DataFrame rename_ped = Rcpp::clone(ped);
	rename_ped[id_name] = int_id;
	rename_ped[dam_name] = int_dam;
	rename_ped[sire_name] = int_sire;
	rename_ped.attr("class") = "data.frame";

	CharacterVector rename_names = as<CharacterVector>(rename_ped.names());
	std::vector<std::string> std_names = rename_ped.names();
	if (id != id_col) rename_names[get_index(std_names, id_name)] = id_col;
	if (dam != dam_col) rename_names[get_index(std_names, dam_name)] = dam_col;
	if (sire != sire_col) rename_names[get_index(std_names, sire_name)] = sire_col;
	rename_ped.attr("names") = rename_names;

	return rename_ped;
}

//' Individuals to be evaluated in purging analyses
//'
//' Returns a boolean vector indicating what individuals are suitable for purging analyses, given a measure of fitness.
//' Individuals with NA values of fitness, and that do not have descendants with non-NA fitness values, are excluded.
//' 
//' @template ped-arg
//' @template valuefrom-arg
//' @return Boolean vector indicating what individuals will be evaluated.
// [[Rcpp::export]]
LogicalVector evaluate(DataFrame ped,
                       std::string value_from) {

  // Read fitness
  Rcpp::CharacterVector w = ped[value_from];

  // Initialize index of IDs to be returned
  R_xlen_t nRows = ped.rows();
  Rcpp::LogicalVector eval (nRows, true);

  // Last individuals with NA data are removed
  // (they do not contribute to inbreeding of any individual with non-NA fitness)
  int last_idx (nRows-1);
  for (int i(last_idx); i>0; --i) {
    if (!Rcpp::CharacterVector::is_na(w(i))) {
      last_idx = i;
      break;
    }
  }
  for (int i(last_idx+1); i<nRows; ++i) eval[i] = false;

  // Remove other individuals with NA fitness and without descendant
  nRows = sum(eval);
  Rcpp::IntegerVector id = ped[id_col];
  Rcpp::IntegerVector dam = ped[dam_col];
  Rcpp::IntegerVector sire = ped[sire_col];
  std::vector<int> dams;
  std::vector<int> sires;
  for (int i(nRows-1); i>0; --i) {
    int id_i = id(i);
    if (Rcpp::CharacterVector::is_na(w(i)) && ((!is_in(id_i, dams)) && (!is_in(id_i, sires)))) eval[i] = false;
    dams.push_back(dam(i));
    sires.push_back(sire(i));
  }
  return eval;
}
