import pandas as pd
from pyproj import Transformer

class DadosEquipamentosPublicos:
    def __init__(self, path):
        
        self.data_frame = pd.read_csv(path, sep=',', encoding='utf-8')
        
        self.data_frame = self.data_frame.drop(columns=['FID', 'link', 'Telefone','DDD','exibirWeb']).reset_index(drop=True)
        
        self.data_frame.rename(columns={
            'geomEquipamentoPublico': 'coordenadas',
        }, inplace=True)
        
        transformer = Transformer.from_crs("EPSG:31983", "EPSG:4326", always_xy=True)
        
        # Função para transformar as coordenadas
        def transform_coordinates(geom):
            # Remove "POINT" e os parênteses, e divide as coordenadas
            x, y = map(float, geom.strip("POINT ()").split())
            # Transforma as coordenadas
            return transformer.transform(x, y)
        
        # Aplica a transformação e cria a tupla
        self.data_frame['coordenadas'] = self.data_frame['coordenadas'].apply(transform_coordinates)
        self.data_frame['coordenadas'] = self.data_frame['coordenadas'].apply(lambda coords: (float(coords[1]), float(coords[0])))

    def get_data(self):
        return self.data_frame
    
    def get_columns(self):
        return self.data_frame.columns.tolist()
    
    def get_schools(self):
        df = self.data_frame[self.data_frame['tematica'] == 'Educação'].reset_index(drop=True)     
            
        return df
    

