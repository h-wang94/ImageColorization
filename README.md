IDEAS:
-Gabor filter parametrizado pelo tamanho da imagem
-Edges no relabeling
-Clustering: Incluir um canal de textura (emular a classificação manual)
-Dimens Reduct: Encontrar melhor subspace
-GP-criptor ? (custo muito elevado)

	
TODO:
-Pasta separada dentro do src para códigos auxiliares
-Color Refinement: Ajuste sat do código de Gupta


Code:
-Colocar a extração de superpixel no final
	-Seguindo todos os ajustes (labeling, sampling, feature extraction)
	-Assim mais fácil de manter o código das outras versões (e os plots intermediários)
-Padronizar as estruturas do código (fluxos de código)
	-source, target, samples, clusters 
-Explicitar os parâmetros escondidos
	-Quando passo a estrutura inteira para uma função e uso internamente os atributos.
