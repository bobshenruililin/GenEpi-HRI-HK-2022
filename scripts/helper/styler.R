# install if you don’t have it
install.packages("styler")

# run this once:
styler::style_file(
  "scripts/2.tree_visualization.R",
  indent_by = 2
)
