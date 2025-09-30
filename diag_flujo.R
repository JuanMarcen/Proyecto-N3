rm(list = ls())
# install.packages("flowchart")
# install.packages('diagram')
# install.packages("DiagrammeR")
# library(flowchart)
# library(diagram)
library(DiagrammeR)
# 
# grViz("
# digraph flujo {
#   graph [layout = dot, rankdir = LR]
#   
#   node [shape = ellipse, style = filled, fillcolor = lightgreen]
#   Inicio
#   
#   node [shape = rectangle, fillcolor = lightblue]
#   Paso1 -> Paso2
#   
#   node [shape = diamond, fillcolor = orange]
#   Decision
#   
#   node [shape = ellipse, fillcolor = lightcoral]
#   Fin
#   
#   Inicio -> Paso1 -> Decision
#   Decision -> Paso2 [label = 'Sí']
#   Decision -> Fin [label = 'No']
#   Paso2 -> Fin
# }
# ")


library(DiagrammeR)


graph <- create_graph()

#1
graph <- graph %>% 
  add_node(
    label = 'Ajuste GLM',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  )

#2
graph <- graph %>% 
  add_node(
    label = 'Validado',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  )

graph <- graph %>%
  add_edge(
    from = 1, to = 2,
    edge_aes = edge_aes(
      color = "black"
    )
  )

#3
graph <- graph %>% 
  add_node(
    label = 'IDK',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  )

graph <- graph %>%
  add_edge(
    from = 2, to = 3,
    edge_aes = edge_aes(
      color = "black",
      label = 'SI'
    )
  )

#4
graph <- graph %>% 
  add_node(
    label = 'IDK',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  )

graph <- graph %>%
  add_edge(
    from = 2, to = 4,
    edge_aes = edge_aes(
      color = "black",
      label = 'NO'
    )
  )

#5
graph <- graph %>%
  add_node(
  label = 'Selección \n covariables',
  type = "type_a",
  node_aes = node_aes(
    shape = 'rectangle',
    color = "steelblue",
    fillcolor = "lightblue",
    fontcolor = "black",
    width = 1,
    height = 0.5
  ),
  node_data = node_data(
    value = 2.5
  )
) %>%
  add_edge(
    from = 1, to = 5,
    edge_aes = edge_aes(
      color = "black"
    )
  )

#6 
graph <- graph %>%
  add_node(
    label = 'glm final',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  ) %>%
  add_edge(
    from = 5, to = 6,
    edge_aes = edge_aes(
      color = "black"
    )
  )%>%
  add_edge(
    from = 6, to = 1,
    edge_aes = edge_aes(
      color = "black"
    )
  )


#7 
graph <- graph %>%
  add_node(
    label = 'OCURRENCIA \n DIARIA',
    type = 'type_a',
    node_aes = node_aes(
      shape = 'rectangle',
      color = "orange",
      fillcolor = "#FCEAD9",
      fontcolor = "black",
      width = 1,
      height = 0.5
    )
  ) %>% 
  add_edge(
    from = 1, to = 7,
    edge_aes = edge_aes(
      color = "black"
    )
  )

#8
graph <- graph %>%
  add_node(
    label = 'Llueve',
    type = 'type_a',
    node_aes = node_aes(
      shape = 'polygon',
      color = "orange",
      fillcolor = "#FCEAD9",
      fontcolor = "black",
      width = 1,
      height = 0.5,
      sides = 4,
      orientation = 45
    )
  ) %>% 
  add_edge(
    from = 7, to = 8,
    edge_aes = edge_aes(
      color = "black"
    )
  )


#9
graph <- graph %>%
  add_node(
    label = 'Lluvia diaria = 0l',
    type = 'type_a',
    node_aes = node_aes(
      shape = 'rectangle',
      color = "orange",
      fillcolor = "#FCEAD9",
      fontcolor = "black",
      width = 2,
      height = 0.5
    )
  ) %>% 
  add_edge(
    from = 8, to = 9,
    edge_aes = edge_aes(
      label = 'No',
      color = "black"
    )
  )

#8
graph <- graph %>%
  add_node(
    label = 'CANTIDAD \n DIARIA',
    type = 'type_a',
    node_aes = node_aes(
      shape = 'rectangle',
      color = "orange",
      fillcolor = "#FCEAD9",
      fontcolor = "black",
      width = 1,
      height = 0.5
    )
  ) %>% 
  add_edge(
    from = 8, to = 10,
    edge_aes = edge_aes(
      label = 'Si',
      color = "black"
    )
  )


#11
graph <- graph %>% 
  add_node(
    label = 'Ajuste GAMLSS',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  ) %>%
  add_edge(
    from = 11, to = 10,
    edge_aes = edge_aes(
      color = "black"
    )
  )

#12
graph <- graph %>% 
  add_node(
    label = 'Validado',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    )
  ) %>%
  add_edge(
    from = 11, to = 12,
    edge_aes = edge_aes(
      color = "black"
    )
  )



#13
graph <- graph %>% 
  add_node(
    label = 'IDK',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  ) %>%
  add_edge(
    from = 12, to = 13,
    edge_aes = edge_aes(
      color = "black",
      label = 'SI'
    )
  )



#14
graph <- graph %>% 
  add_node(
    label = 'IDK',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    ),
    node_data = node_data(
      value = 2.5
    )
  ) %>%
  add_edge(
    from = 12, to = 14,
    edge_aes = edge_aes(
      color = "black",
      label = 'NO'
    )
  )

#15
graph <- graph %>%
  add_node(
    label = 'Selección \n covariables',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    )
  ) %>%
  add_edge(
    from = 11, to = 15,
    edge_aes = edge_aes(
      color = "black"
    )
  )

#16 
graph <- graph %>%
  add_node(
    label = 'GAMLSS final',
    type = "type_a",
    node_aes = node_aes(
      shape = 'rectangle',
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "black",
      width = 1,
      height = 0.5
    )
  ) %>%
  add_edge(
    from = 15, to = 16,
    edge_aes = edge_aes(
      color = "black"
    )
  )%>%
  add_edge(
    from = 16, to = 11,
    edge_aes = edge_aes(
      color = "black"
    )
  )


render_graph(graph, layout = 'nicely')

# grViz("
# digraph G {
#   
#   # Layout
#   graph [rankdir = TB]
# 
#   # Nodos principales
#   AjusteGLM [shape=rectangle, style=filled, fillcolor=lightblue, label='Ajuste GLM']
#   Validado [shape=rectangle, style=filled, fillcolor=lightblue, label='Validado']
# 
#   # Flechas principales
#   AjusteGLM -> Validado
#   AjusteGLM -> MDOCluster
# 
#   # Nodo grande como cluster
#   subgraph cluster_MDO {
#     label = 'MDO';         # título del cluster
#     style = rounded;        # borde redondeado
#     color = steelblue;      # color del borde del cluster
#     node [style=filled, fillcolor=lightyellow]
# 
#     MDO_1 [label='Nodo 1']
#     MDO_2 [label='Nodo 2']
# 
#     MDO_1 -> MDO_2
#   }
# 
# }
")
