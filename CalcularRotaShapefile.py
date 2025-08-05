# Importa bibliotecas necessárias para geoprocessamento e construção de grafos
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
from pyproj import Transformer
from shapely.strtree import STRtree
import numpy as np
from scipy.spatial import KDTree
import random
from shapely.ops import split

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
        
        self.gdf = self.gdf.explode().reset_index(drop=True).copy()
        
        self.gdf["grafo"] = self.gdf["geometry"].apply(lambda geometry: [self._arredondar(x) for x in geometry.coords])
        self._inserir_nos_em_intersecoes()
        # Constrói o grafo com base nas linhas da rede viária
        self.G = self._criar_grafo()

        # Define transformador para converter coordenadas UTM → WGS84 (lat/lon)
        self._transformador_para_leaflet = Transformer.from_crs(self.gdf.crs, "EPSG:4326", always_xy=True)
        
        
    def _arredondar(self,coord, precisao=10):
        """
        Arredonda uma coordenada para reduzir duplicidade de nós.
        """
    
        
        return (round(coord[0], precisao), round(coord[1], precisao))

    def _criar_grafo(self):
        """
        Constrói o grafo a partir das linhas do GeoDataFrame.
        """
        G = nx.Graph()  # Cria grafo não-direcionado
        for index, linha in self.gdf.iterrows():
            # Extrai as coordenadas das linhas
            coords = list(linha.geometry.coords)
            # Adiciona arestas ao grafo com pesos baseados na distância
        
            for i in range(len(coords) - 1):

                u = linha["grafo"][i]
                v = linha["grafo"][i+1]
                distancia = Point(coords[i]).distance(Point(coords[i + 1]))
                G.add_edge(u, v, distancia=distancia)
       
        #G = self._unir_componentes_desconectados(G)
        return G
    

    def _inserir_nos_em_intersecoes(self):
        """
        Insere nós nas interseções entre as ruas, dividindo as LineStrings
        para garantir que os cruzamentos virem vértices no grafo.
        """
        linhas = list(self.gdf.geometry)
        arvore = STRtree(linhas)
        
        novas_linhas = []
        for linha in linhas:
            # Coleta todas as geometrias que cruzam a linha atual
            intersecoes = []
            for outra in arvore.query(linha):
                if linha == outra:
                    continue
                print(f'outra: {outra} linha: {linha}')
                if linha.crosses(outra) or linha.touches(outra) or linha.intersects(outra):
                    ponto = linha.intersection(outra)
                    if ponto.geom_type == 'Point':
                        intersecoes.append(ponto)
                    elif ponto.geom_type == 'MultiPoint':
                        intersecoes.extend(ponto.geoms)

            # Divide a linha nas interseções
            if intersecoes:
                for ponto in intersecoes:
                    if ponto.within(linha):
                        linha = split(linha, ponto)
                        # Se o split retorna uma coleção de linhas, seguimos
                        if hasattr(linha, 'geoms'):
                            linha = list(linha.geoms)
                        else:
                            linha = [linha]
                # Se foi dividida, adiciona cada segmento
                novas_linhas.extend(linha)
            else:
                novas_linhas.append(linha)

        # Substitui o GeoDataFrame com as novas linhas segmentadas
        self.gdf = gpd.GeoDataFrame(geometry=novas_linhas, crs=self.gdf.crs)

    def _unir_componentes_desconectados(self, G):
        """
        Une componentes desconectados de um grafo ligando os nós mais próximos entre componentes.
        O grafo deve ter os nós como tuplas de coordenadas (x, y).
        """

        while nx.number_connected_components(G) > 1:
            componentes = list(nx.connected_components(G))
            print(nx.number_connected_components(G))
            # Cria um KDTree para cada componente e guarda os nós
            kdtrees = []
            coord_nos = []
            for comp in componentes:
                nos = list(comp)
                coords = np.array(nos)  # Os nós são as próprias coordenadas
                kdtrees.append(KDTree(coords))
                coord_nos.append(nos)

            # Busca o menor par de nós entre qualquer par de componentes
            menor_dist = float("inf")
            melhor_par = (None, None)

            for i in range(len(componentes)):
                for j in range(i + 1, len(componentes)):
                    dists, idxs = kdtrees[i].query(np.array(coord_nos[j]))
                    min_idx = np.argmin(dists)
                    dist = dists[min_idx]
                    if dist < menor_dist:
                        menor_dist = dist
                        no1 = coord_nos[i][idxs[min_idx]]
                        no2 = coord_nos[j][min_idx]
                        melhor_par = (no1, no2)

            # Conecta os dois nós mais próximos entre componentes diferentes
            if melhor_par[0] and melhor_par[1]:
                distancia = Point(melhor_par[0]).distance(Point(melhor_par[1]))
                print(f"Distancia da juntação: {distancia}")
                G.add_edge(melhor_par[0], melhor_par[1], distancia=distancia)

        return G

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

    def gerar_componentes_com_cores_para_leaflet(self):
        """
        Gera uma lista de componentes do grafo, cada um com uma cor distinta e suas coordenadas convertidas
        para latitude/longitude, para visualização com Leaflet.
        
        :return: Lista de dicionários com chaves: 'cor' e 'coordenadas' (lista de pares lat/lon).
        """
        componentes = list(nx.connected_components(self.G))
        resultados = []

        for componente in componentes:
            subgrafo = self.G.subgraph(componente)
            coordenadas_componentes = []

            for u, v in subgrafo.edges():
                # Cada aresta é uma linha entre dois pontos
                lon1, lat1 = self._transformador_para_leaflet.transform(u[0], u[1])
                lon2, lat2 = self._transformador_para_leaflet.transform(v[0], v[1])
                coordenadas_componentes.append([(lat1, lon1), (lat2, lon2)])

            # Gera uma cor aleatória em formato hexadecimal
            cor_hex = "#{:06x}".format(random.randint(0, 0xFFFFFF))

            resultados.append({
                "cor": cor_hex,
                "coordenadas": coordenadas_componentes
            })

        return resultados