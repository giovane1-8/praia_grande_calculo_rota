# Importa bibliotecas necessárias para geoprocessamento e construção de grafos
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
from pyproj import Transformer
from shapely.strtree import STRtree

class CalcularRotaShapefile:
    def __init__(self, caminho_shp: str):
        """
        Inicializa a classe com base em um shapefile da rede viária.
        :param caminho_shp: caminho para o shapefile com a rede de caminhos (linhas) em EPSG:31983.
        """
        self.gdf = gpd.read_file(caminho_shp)                            # Lê o shapefile como GeoDataFrame
        
        # Verifica se o sistema de coordenadas está correto (UTM zona 23S)
        if self.gdf.crs.to_epsg() != 31983:
            raise ValueError("O shapefile deve estar no sistema de coordenadas EPSG:31983 (UTM 23S)")
        
        # Remove vias não apropriadas para pedestres (ex: rodovias e vias expressas)
        self.gdf = self.gdf[~self.gdf["codigoTipo"].isin([6, 8])].copy() 
        
        self.gdf = self.gdf.explode().reset_index(drop=True)
        

        # Constrói o grafo com base nas linhas da rede viária
        self.G = self._criar_grafo()

        # Define transformador para converter coordenadas UTM → WGS84 (lat/lon)
        self._transformador_para_leaflet = Transformer.from_crs(self.gdf.crs, "EPSG:4326", always_xy=True)

    def _criar_grafo(self):
        """
        Constrói o grafo a partir das linhas do GeoDataFrame.
        """
        G = nx.Graph()  # Cria grafo não-direcionado

        for geom in self.gdf.geometry:
            # Extrai as coordenadas das linhas
            coords = list(geom.coords)
            # Adiciona arestas ao grafo com pesos baseados na distância
            for i in range(len(coords) - 1):

                u = self._arredondar(coords[i])
                v = self._arredondar(coords[i + 1])
                distancia = Point(coords[i]).distance(Point(coords[i + 1]))
                G.add_edge(u, v, distancia=distancia)

        return G
  
        
        
    def _arredondar(self, coord, precisao=0):
        """
        Arredonda uma coordenada para reduzir duplicidade de nós.
        """
        return (round(coord[0], precisao), round(coord[1], precisao))

    def _inserir_no_temporario(self, coord_utm, grafo):
        """
        Insere um nó temporário projetado sobre a aresta mais próxima do grafo.
        Se não encontrar aresta, adiciona como nó isolado.
        """
        ponto_proj = self._arredondar(coord_utm)
        ponto_geom = Point(ponto_proj)

        # Cria lista de geometrias das arestas com seus dados
        arestas_com_geometrias = []
        for u, v, data in grafo.edges(data=True):
         
            if not isinstance(u, tuple) or not isinstance(v, tuple):
                continue
            if u == v:
                continue
            arestas_com_geometrias.append((LineString([u, v]), (u, v, data)))

        if not arestas_com_geometrias:
            # Caso não haja arestas válidas
            grafo.add_node(ponto_proj, x=ponto_proj[0], y=ponto_proj[1])
            return ponto_proj

        # Cria índice espacial para busca rápida de linhas próximas
        tree = STRtree([line_obj for line_obj, _ in arestas_com_geometrias])

        menor_distancia = float("inf")
        aresta_mais_proxima = None
        ponto_projetado_mais_proximo = None
        candidate_arestas = []

        # Encontra a linha mais próxima usando o índice espacial
        try:
            closest_line_geom_indices = tree.query_nearest(ponto_geom, all_matches=False)

            if closest_line_geom_indices.size > 0:
                closest_line_geom_idx_scalar = closest_line_geom_indices[0]
                closest_line_geom, original_edge_data = arestas_com_geometrias[closest_line_geom_idx_scalar]
                candidate_arestas.append((closest_line_geom, original_edge_data))
        except AttributeError:
            # Se não suportar query_nearest (ex: Shapely < 2.0)
            candidate_arestas = arestas_com_geometrias

        # Itera nas arestas candidatas e encontra a projeção mais próxima
        for linha, (u, v, data) in candidate_arestas:
            projecao = linha.interpolate(linha.project(ponto_geom))
            projecao_tuple = self._arredondar((projecao.x, projecao.y))

            dist_total = linha.length
            dist_u = Point(u).distance(projecao)
            dist_v = Point(v).distance(projecao)
            soma = dist_u + dist_v

            # Verifica se o ponto projetado realmente pertence ao segmento
            if abs(soma - dist_total) < 1e-3:
                distancia_ao_segmento = ponto_geom.distance(projecao)
                if distancia_ao_segmento < menor_distancia:
                    menor_distancia = distancia_ao_segmento
                    aresta_mais_proxima = (u, v)
                    ponto_projetado_mais_proximo = projecao_tuple

        # Se encontrou uma aresta próxima válida
        if aresta_mais_proxima:
            u, v = aresta_mais_proxima

            if grafo.has_edge(u, v):
                grafo.remove_edge(u, v)  # Remove aresta original

            # Adiciona novo nó projetado
            grafo.add_node(ponto_projetado_mais_proximo,
                           x=ponto_projetado_mais_proximo[0],
                           y=ponto_projetado_mais_proximo[1])

            # Conecta novo nó aos nós originais com as distâncias corretas
            dist_u = Point(u).distance(Point(ponto_projetado_mais_proximo))
            dist_v = Point(v).distance(Point(ponto_projetado_mais_proximo))
            if ponto_projetado_mais_proximo != u:
                grafo.add_edge(u, ponto_projetado_mais_proximo, distancia=dist_u)
            if ponto_projetado_mais_proximo != v:
                grafo.add_edge(ponto_projetado_mais_proximo, v, distancia=dist_v)

            return ponto_projetado_mais_proximo

        # Se não encontrou nenhuma aresta válida, adiciona como nó isolado
        grafo.add_node(ponto_proj, x=ponto_proj[0], y=ponto_proj[1])
        return ponto_proj

    def calcular_rota(self, origem: tuple, destino: tuple):
        """
        Calcula a menor rota entre dois pontos geográficos (lat, lon).
        :param origem: tupla (lat, lon)
        :param destino: tupla (lat, lon)
        :return: (distancia_em_metros, lista de pontos [(lat, lon), ...])
        """
        grafo = self.G.copy()  # Cria uma cópia do grafo original

        # Converte coordenadas WGS84 → UTM
        transformador_para_31983 = Transformer.from_crs("EPSG:4326", self.gdf.crs, always_xy=True)
        origem_proj = transformador_para_31983.transform(origem[1], origem[0])
        destino_proj = transformador_para_31983.transform(destino[1], destino[0])

        # Insere nós temporários na origem e destino
        no_origem = self._inserir_no_temporario(origem_proj, grafo)
        no_destino = self._inserir_no_temporario(destino_proj, grafo)

        # Calcula o menor caminho usando o peso 'weight' (distância)
        try:    
            path = nx.shortest_path(grafo, no_origem, no_destino, weight="distancia")
            dist = nx.path_weight(grafo, path, weight="distancia")
        except:
            # Se não houver caminho possível
            path = [no_destino]
            dist = -9999

        rota = []
        for x, y in path:
            # Converte coordenadas para (lat, lon) para exibir no mapa
            lon, lat = self._transformador_para_leaflet.transform(x, y)
            rota.append((lat, lon))

        return dist, rota  # Retorna a distância e o caminho
