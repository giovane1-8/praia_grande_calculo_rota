from flask import Flask, request, render_template, jsonify

from CalcularRotaShapefile import CalcularRotaShapefile
from DadosEquipamentosPublicos import DadosEquipamentosPublicos 

 
app = Flask(__name__, template_folder='./')


equipamentos_publicos = DadosEquipamentosPublicos('EquipamentoPublicoDados.csv')


calcular_rota = CalcularRotaShapefile("MapaTrechoCidade/MapaTrechoCidadeLine.shp")

@app.route('/', methods=['GET', 'POST'])
def sua_rota_de_salvar():
    if request.method == 'POST':
        latitude = float(request.form['latitude'])
        longitude = float(request.form['longitude'])
        resultado = []
        escolas = equipamentos_publicos.get_schools()
       
        for index, escola in escolas.iterrows():
        
            origem = (latitude,longitude)
           
            destino = escola['coordenadas']
         
            print(f"{index}: {escola['Nome']}")
            distancia, rota = calcular_rota.calcular_rota(origem, destino)  # rota deve ser lista de [lat, lon]
            resultado.append({
                "escola": escola['Nome'],
                "distancia": distancia,
                "rota": rota,
                "latitude": escola['coordenadas'][0],
                "longitude": escola['coordenadas'][1]
            })

        resultado.sort(key=lambda x: x['distancia'])


        return render_template('resultado.html', resultado=resultado, origem=( latitude, longitude))

    return render_template('cadastro.html')
@app.route('/api/<origem>/<destino>/', methods=['GET'])
def api_calcular_rota(origem, destino):
    # Converte as coordenadas de origem e destino para tuplas
    try: 
        origem = tuple(map(float, origem.split(',')))
        destino = tuple(map(float, destino.split(',')))
    except ValueError:
        return jsonify({"error": "Coordenadas inválidas"}), 400
    # Chama o método de cálculo de rota
    distancia, rota = calcular_rota.calcular_rota(origem, destino)

    # Retorna a resposta em formato JSON
    return jsonify({
        "origem": origem,
        "destino": destino,
        "distancia": distancia,
        "rota": rota
    })



@app.route('/verificarcomponentes/', methods=['GET'])
def mapa():
    return render_template('mostra_grafo_mapa.html', componentes=calcular_rota.gerar_componentes_com_cores_para_leaflet())