library(DiagrammeR)

#----datos---- 
# graph_datos <- create_graph()
# graph_datos <- graph_datos %>%
#   add_node(
#     label = "Datos 15min",
#     node_aes = node_aes(
#       shape = 'rectangle',
#       color = "steelblue",
#       fontcolor = "black",
#       width = 1,
#       height = 0.5,
#       x = 0,
#       y = 0
#     )
#   ) %>%
#   add_node(
#     label = 'Datos horarios',
#     node_aes = node_aes(
#       shape = 'rectangle',
#       color = "steelblue",
#       fontcolor = "black",
#       width = 1.5,
#       height = 0.5,
#       x = 2.5,
#       y = 0
#     )
#   ) %>%
#   add_node(
#     label = 'Datos diarios',
#     node_aes = node_aes(
#       shape = 'rectangle',
#       color = "steelblue",
#       fontcolor = "black",
#       width = 1.5,
#       height = 0.5,
#       x = 5,
#       y = 0
#     )
#   ) %>%
#   add_edge(
#     from = 1, to = 2,
#     edge_aes = edge_aes(
#       color = "black", 
#       label = "15min NA --> NA"
#     )
#   ) %>%
#   add_edge(
#     from = 2, to = 3,
#     edge_aes = edge_aes(
#       color = "black",
#       label = '+1/4h NA --> NA'
#     )
#   )
# render_graph(graph_datos, layout = 'tree')


#----generador diario----

gd <- create_graph()
gd <- gd %>%
  add_node(
    label = "Datos diarios",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  )%>%
  add_node(
    label = "Datos con \n existencias",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_edge(
    from = 1, to = 2,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

gd <- gd %>%
  add_node(
    label = "Tratamiento \n GLM (logistic reg.)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>% 
  add_node(
    label = "Lag lluvia \n diaria",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = " Variables climáticas \n G, T",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fontcolor = "black",
      width = 1.5,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Armónicos",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  )

gd <- gd %>%
  add_edge(
    from = 2, to = 3,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 4, to = 3,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 5, to = 3,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 6, to = 3,
    edge_aes = edge_aes(
      color = 'black'
    )
  )


gd <- gd %>%
  add_node(
    label = "Estudio transf. \n poly variables",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Selección AIC \n variables",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_edge(
    from = 3, to = 7,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 7, to = 8,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

gd <- gd %>% 
  add_node(
    label = "GLM final",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>% 
  add_node(
    label = "Ajuste GLM",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>% 
  add_edge(
    from = 8, to = 9,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 9, to = 10,
    edge_aes = edge_aes(
      color = 'black'
    )
  )
  
gd <- gd %>%
  add_node(
    label = "Validación",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "??",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "??",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "blue",
      fillcolor = "#BAC1FF",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>% 
  add_edge(
    from = 10, to = 11,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>% 
  add_edge(
    from = 11, to = 12,
    edge_aes = edge_aes(
      color = 'black', 
      label = 'SÍ'
    )
  ) %>% 
  add_edge(
    from = 11, to = 13,
    edge_aes = edge_aes(
      color = 'black', 
      label = 'NO'
    )
  )

gd <- gd  %>%
  add_node(
    label = "Ocurrencia Diaria \n (MDO)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  )  %>%
  add_edge(
    from = 10, to = 14,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

gd <- gd %>%
  add_node(
    label = "FOR day 1 to end",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Muestra U(0,1)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "U(0,1) ≤ P(llueve.day = 1 | datos)",
    node_aes = node_aes(
      shape = 'polygon',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 2,
      height = 1,
      sides = 4,
      orientation = 45
    )
  ) %>%
  add_edge(
    from = 14, to = 15,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 15, to = 16,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 16, to = 17,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

gd <- gd  %>%
  add_node(
    label = "??",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Cantidad diaria \n (MDQ)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_edge(
    from = 17, to = 18,
    edge_aes = edge_aes(
      color = 'black',
      label = 'NO'
    )
  ) %>%
  add_edge(
    from = 17, to = 19,
    edge_aes = edge_aes(
      color = 'black',
      label = 'SI'
    )
  )

gd <- gd %>%
  add_edge(
    from = 13, to = 18,
    edge_aes = edge_aes(
      color = 'black',
      style = 'dashed'
    )
  )

gd <- gd  %>%
  add_node(
    label = "Predict shape r and rate λ \n of Gamma(r, λ)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.5,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "LLuvia.day = sample Gamma(r, λ)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#C72100",
      fillcolor = "#F29380",
      fontcolor = "black",
      width = 1.5,
      height = 0.5,
    )
  ) %>%
  add_edge(
    from = 19, to = 20,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 20, to = 21,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

#ajuste MDQ
gd <- gd %>%
  add_node(
    label = "Tratamiento \n  GAMLSS(glm gamma)",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Estudio trans \n poly",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Selección AIC \n variables",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  )  %>%
  add_node(
    label = "GAMLSS final",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Ajuste GAMLSS",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Validación",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "??",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "??",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  ) %>%
  add_node(
    label = "Datos lluvia > 0",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "#009E00",
      fillcolor = "#C9FFC9",
      fontcolor = "black",
      width = 1.25,
      height = 0.5,
    )
  )


gd <- gd %>%
  add_edge(
    from = 22, to = 23,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 23, to = 24,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 24, to = 25,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 25, to = 26,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 26, to = 27,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 27, to = 28,
    edge_aes = edge_aes(
      color = 'black',
      label = 'SÍ',
    )
  ) %>%
  add_edge(
    from = 27, to = 29,
    edge_aes = edge_aes(
      color = 'black', 
      label = 'NO'
    )
  ) %>%
  add_edge(
    from = 26, to = 19,
    edge_aes = edge_aes(
      color = 'black'
    )
  )

#las de tratamiento
gd <- gd  %>%
  add_edge(
    from = 4, to = 22,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 5, to = 22,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 6, to = 22,
    edge_aes = edge_aes(
      color = 'black'
    )
  ) %>%
  add_edge(
    from = 30, to = 22,
    edge_aes = edge_aes(
      color = 'black'
    )
  )


#vuelta bucle 
gd <- gd %>%
  add_edge(
    from = 21, to = 15,
    edge_aes = edge_aes(
      color = 'black'
    )
  )




render_graph(gd)
export_graph(graph = gd, file_name = "mi_grafico.svg", file_type = "SVG")

