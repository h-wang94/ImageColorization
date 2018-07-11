IDEAS:
-Clustering de features no target (spatial coherence)
	-Mesmo cluster + conectado = mesma classe. (voto majoritário)
-Gabor filter parametrizado pelo tamanho da imagem
-Superpixels: Feature Averaging (média pode ter outliers)
	-Ao invés de averaging, guardar vetores.
	-Distância entre superpixels definida como a média entre as menores distâncias em pixel
-Edges no relabeling
-Clustering: Incluir um canal de textura (emular a classificação manual)
-Dimens Reduct: Encontrar melhor subspace
-GP-criptor ? (custo muito elevado)

	
TODO:
-Pasta separada dentro do src para código auxiliar
-Relabeling: Ajuste dos parâmetros
-Color Refinement: Ajuste sat do código de Gupta
-Corrigir todas as chamadas à função "mode"
	-especialmente a do label automático
-Padronizar nomes dos atributos do source e target.
-Features restantes


Code:
-Criar uma struct para os superpixels (análoga a struct 'samples')
*Repensar a estrutura (source, target, samples e superpixels)

