# Function to draw a valuebox in bare Shiny (maybe replace with dedicated dashboard packages later)

valueBox <- function(value, label, color = "#0073b7", width = "200px") {
  div(
    class = "value-box",
    style = paste0(
      "display:inline-block; ",
      "padding:20px; margin:10px; border-radius:8px; ",
      "background-color:", color, "; color:white; ",
      "width:", width, "; text-align:center; font-weight:bold; ",
      "box-shadow:2px 2px 6px rgba(0,0,0,0.2);"
    ),
    div(style = "font-size:32px;", value),
    div(style = "font-size:16px; opacity:0.8;", label)
  )
}