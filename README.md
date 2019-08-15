# optimization-algorithms-on-matrix-manifolds-with-R
用R语言将&lt;optimization algorithms on matrix manifolds>这本书内的几个算法实现了一下

+ algorithm3: 线性搜索算法
+ algorithm8: 牛顿法
+ algorithm10&12: 信赖域法
+ algorithm13: 共轭梯度法

解决的都是Grass(p,n)上的瑞利商的最优化，也就是矩阵的特征值与广义特征值问题

注：algorithm8貌似只能收敛到局部最优点

I have apply some algorithms of the book <optimization algorithms on matrix manifold> with R 

+ algorithm3: line search method
+ algorithm8: Newton method
+ algorithm10&12: trust region method
+ algorithm13: conjucated gradient method

All of the algorithms are used to solve rayleigh quotient problem on Grass(p, n), which can also be transformed to 
eigenvalue or generalized eigenvalue problems

ps: it seems that algorithm8 can only converges to local convergence points
