linters <- lintr::linters_with_defaults(
  line_length_linter = NULL,
  commented_code_linter = NULL,
  trailing_whitespace_linter = NULL,
  object_name_linter = lintr::object_name_linter(styles = "camelCase"),
  return_linter = lintr::return_linter(return_style = "explicit"),
  seq_linter = NULL
)
