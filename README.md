IDEAS:
-Gabor filter parametrizado pelo tamanho da imagem
-nClusters (cor) ajustado pelo erro de recolorização.
-Edges no relabeling
-Clustering: Incluir um canal de textura (emular a classificação manual)
-GP-criptor ? (custo muito elevado)

	
TODO:
-Refactoring:
	-Parei em: Source sampling
-Pasta separada dentro do src para códigos auxiliares
-Color Refinement: Ajuste sat do código de Gupta


Code:
-Dimensionality Reduction por feature
	-Criar vetor de ativação de Dim Reduction
	-Incluir no batch da busca (?)
-Tirar os parâmetros Hard (especialmente o valor 40 do peaks).
-Colocar a extração de superpixel no final
	-Seguindo todos os ajustes (labeling, sampling, feature extraction)
	-Assim mais fácil de manter o código das outras versões (e os plots intermediários)
-Padronizar as estruturas do código (fluxos de código)
	-source, target, samples, clusters 
-Explicitar os parâmetros escondidos
	-Quando passo a estrutura inteira para uma função e uso internamente os atributos.
-Indexar todas as chamadas de "figure"