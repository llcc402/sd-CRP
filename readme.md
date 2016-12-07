# The Spectral Dependent Chinese Restaurant Process

This method is proposed by Socher, Richard; Maas, Andrew; Manning, Christopher D. Published by Proceedings of the Fourteenth International Conference on Artificial Intelligence and Statistics, 2011(15), p 698 - 706.

## Discription of sd-CRP

It is a generalization of Blei's Distance Dependent Chinese Restaurant Process.

It calculates the distance between every pair of observations using the Laplacian transformation. Concretely, let L denotes the NORMALIZED Laplacian matrix of the data set, U is a matrix composed of the eigenvectors of L, whose order is m by n, where m is the number of the observations and n is a presumed dimension of U. Then the distance of observation i and j is the Euclidean distance between U(i,:) and U(j,:).

With the above distance in hand, it assumes that the probability of observation i to sit with observation j is proportional to

	p(c_i = j) p(x_{1:N}|c_i = j, c_{-i})

where 

	p(c_i = j) = f(d_{ij})

## Functions included

1. data_generate.m:	Generate two classes of observations.

2. get_similarity.m:	Compute the similarity of every pair of obbservations.

3. get_sitBehind.m:	For each observation i, it returns the customers who sitwith i.

4. sdCRP.m:		Returns the table label for every customer.
