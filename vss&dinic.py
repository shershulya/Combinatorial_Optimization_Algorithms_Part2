# x,y-VSS (Vertex separating set или вершинно-разделяющее множество):
# Такое множество вершин что удалив их из исходного графа и инциндентные им рёбра, x и y перестают быть связанными.
# Задача: найти x,y-VSS минимальной мощности
# 
# Строим вспомогательную сеть:
# 1) G = (V', D', x, y)
# 2) Любую вершину v из V\{x,y} заменим двумя вершинами (v, 0) (v, 1)
# 3) Если дуга (v, v') была в исходном графе и v, v'∉ {x, y}, то заменим её на две дуги: ((v, 1)->(v', 0)) и ((v', 1)->(v, 0))
# 4) Дуги (x, v) заменим на (x->(v, 0)) и (v, y) на ((v, 1)->y)
# 
# Если в такой сети найти максимальный поток, то мощность этого максимального потока будет соотвествовать количеству вершин минимального x,y-VSS
# 
# ЭК:
#  | V' | = 2n - 2 = O(n)
#  | D' | = 2m + n = O(m)
#  => Трудоёмкость алгоритма Эдмондса-Карпа: O(m * n^2)
# Диниц:
#  Вершинная характеристика: O(n)
#  Рёберная характеристика: большая, т.к. у каждой вершины есть дуги с бесконечной пропускной способностью
#  => Характеристика: O(n)
#  => Трудоёмкость алгоритма Диница: O((m + n) * √(n)) = O(m * √(n))

import networkx as nx
import matplotlib.pyplot as plt
from copy import deepcopy

class Edge:
  def __init__(self, to, cap, flow, link, color):
    self.to = to
    self.cap = cap
    self.flow = flow
    self.link = link
    self.color = color
    

class Graph:
  def __init__(self, N, S, T):
    self.N = N
    self.S = S
    self.T = T
    self.level = [0 for _ in range(N)]
    self.queue = [0 for _ in range(N)]
    self.head = [-1 for _ in range(N)]
    self.start = [-1 for _ in range(N)]
    self.edges = []
    self.edgelist = []
    
  def AddEdge(self, a, b, cap, color='black'):
    e1 = Edge(b, cap, 0, self.head[a], color)
    e2 = Edge(a, 0, 0, self.head[b], color)
    self.head[a] = len(self.edges)
    self.edges.append(e1)
    self.head[b] = len(self.edges)
    self.edges.append(e2)

    self.edgelist.append((a, b))
    
  def Print(self):
    G = nx.DiGraph()
    G.add_edges_from(self.edgelist)

    layout = {}
    layout[self.S] = (0, 0)
    height = self.T / 4
    for i in range(1, self.N, 2):
      layout[i] = (1, height - i / 2)
      layout[i + 1] = (2, height - i / 2)
    layout[self.T] = (3, 0)
    
    nx.draw_networkx_nodes(G,
      pos = layout, 
      edgecolors='black', 
      linewidths=1, 
      alpha=0.5)

    nx.draw_networkx_labels(G,
      pos=layout,
      labels={i: i for i in range(self.N)})

    nx.draw_networkx_edges(G, 
      pos = layout, 
      alpha=0.5)

    edgelabels_b = {}
    edgelabels_g = {}
    for node in G.nodes():
      id = self.head[node]
      while id != -1:
        e = self.edges[id]
        if e.cap > 0:
          if e.color == 'black':
            edgelabels_b[(node, e.to)] = str(e.flow) + '/' + str(e.cap)
          else:
            edgelabels_g[(node, e.to)] = str(e.flow) + '/' + str(e.cap)
        id = e.link

    nx.draw_networkx_edge_labels(G, 
      pos=layout, 
      edge_labels=edgelabels_b,
      font_color='black',
      label_pos=0.8)

    nx.draw_networkx_edge_labels(G, 
      pos=layout, 
      edge_labels=edgelabels_g,
      font_color='green',
      label_pos=0.8)
    
    plt.show()
    return


  # Finds if more flow can be sent from s to t
  # Also assigns levels to nodes
  def BFS(self):
    self.level = [-1 for _ in range(self.N)]
    self.level[self.S] = 0
    
    # Create a queue, enqueue source vertex
    # and mark source vertex as visited here
    q_curr = 0
    q_new = 1
    self.queue[0] = self.S
    while (q_curr < q_new and self.level[self.T] == -1):
      v = self.queue[q_curr]
      q_curr += 1
      id = self.head[v]
      while id != -1:
        e = self.edges[id]
        to = e.to
        id = e.link
        if self.level[to] == -1 and e.flow < e.cap:
          self.level[to] = self.level[v]+1
          self.queue[q_new] = to
          q_new += 1
          
    # If we can not reach to the sink we return False and Dinic is completed
    return self.level[self.T] != -1
  
  # A DFS based function to send flow after BFS has
  # figured out that there is a possible flow and
  # constructed levels. This functions called multiple
  # times for a single call of BFS.
  def DFS(self, v, flow):
    if not flow:
      return 0
    # Sink reached
    if v == self.T:
      return flow

    id = self.start[v]
    while id != -1:
      # Pick next edge from sheaf list
      e = self.edges[id]
      self.start[v] = e.link
      to = e.to
      if self.level[to] != self.level[v] + 1:
        id = e.link
        continue
      pushed = self.DFS(to, min(flow, e.cap - e.flow))
      if pushed:
        self.edges[id].flow += pushed
        self.edges[id ^ 1].flow -= pushed
        return pushed
      id = e.link
    return 0

  def Dinic(self):
    if self.S == self.T:
      return -1
    flow = 0
    # Augument the flow while there is path
    # from source to sink
    while self.BFS():
      self.start = deepcopy(self.head)
      while pushed := self.DFS(self.S, float('inf')):
        # Add path flow to overall flow
        flow += pushed
    return flow

def VSS(n, s, t, filename):
  g = Graph(2 * n - 2, s, 2 * t - 1)
  with open(filename + ".txt") as file:
    lines = file.readlines()
    for line in lines:
      a, b, cap = map(int, line.split())
      if a == s and b == t or a == t and b == s:
        print("There is no such vertex separating set")
        exit(0)
      if a != s and a != t and b != s and b != t:
        g.AddEdge(2 * a, 2 * b - 1, float('inf'))
        g.AddEdge(2 * b, 2 * a - 1, float('inf'))
      else:
        if a == s:
          g.AddEdge(a, 2 * b - 1, float('inf'))
        if a == t:
          g.AddEdge(2 * b, 2 * a - 1, float('inf'))
        if b == s:
          g.AddEdge(b, 2 * a - 1, float('inf'))
        if b == t:
          g.AddEdge(2 * a, 2 * b - 1, float('inf'))

  for i in range(1, t):
    g.AddEdge(2 * i - 1, 2 * i, 1, 'green')

  print("Maximum flow", g.Dinic())
  g.Print()

if __name__ == "__main__":
  VSS(n = 7, s = 0, t = 6, filename = "test_flow1")
  VSS(n = 6, s = 0, t = 5, filename = "test_flow2")