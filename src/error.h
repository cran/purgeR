#ifndef ERROR_H
#define ERROR_H

#include <stdexcept>
#include <Rcpp.h>

inline void print_error (std::string message, std::string value) {
  throw std::logic_error(message + ": " + value);
}

inline void print_warning (std::string message, std::string value) {
  std::string warning = "Warning! " + message + ": " + value;
  ::Rf_warning(warning.c_str());
}

const std::string

// Unexpected error
ERROR_UNEXPECTED = "Unexpected error, please contact the developer",

// Errors reading options from the terminal
ERROR_INPUT_COL_ID = "Mandatory column 'id' not found",
ERROR_INPUT_COL_DAM = "Mandatory column 'dam' not found",
ERROR_INPUT_COL_SIRE = "Mandatory column 'sire' not found",
ERROR_INPUT_COL_REF_ZERO = "At least one individual must belong to the reference population!",
ERROR_INPUT_IND_ZERO = "Individual IDs cannot be set to 0 (it is reserved to unknown parents)",
ERROR_INPUT_IND_REPEAT = "There are repeated individual IDs",
ERROR_INPUT_IND_CIRCLE = "Individuals cannot be born from themselves!",
ERROR_INPUT_IND_SELFING = "Selfing detected",
ERROR_INPUT_IND_ORDER = "Dams and sires must be declared before their offspring!",
ERROR_INPUT_IND_INDEX = "Individuals must be named from 1 to N",
ERROR_ARGUMENT_ANC_NBOOT = "Number of iterations for bootstrapping must be between 10e2 and 10e8",
ERROR_ARGUMENT_FIJ_ANCESTORS = "Vector of ancestors should only be given with argument mode='custom'",
ERROR_ARGUMENT_FIJ_CUSTOM = "mode='custom' requires a vector of ancestors ('ancestors' argument)",
ERROR_ARGUMENT_FIJ_MODE = "Unknown mode value. Select one of 'founders', 'all', 'custom'",
ERROR_ARGUMENT_INB_D = "Purging coefficient must be in the range [0.0, 0.5]",
ERROR_ARGUMENT_OOE_DEPTH = "Depth must be an integer higher than 0";

// Errors that need to be declared of type char*
//const char* ERROR_INPUT_COL_TYPE = "Mandatory 'id', 'dam' and 'sire' columns need to be of type integer";
const char* ERROR_INPUT_COL_TYPE_W = "Fitness needs to be of numeric type";
const char* ERROR_INPUT_COL_TYPE_R = "Reference needs to be coercible to boolean type";

#endif
