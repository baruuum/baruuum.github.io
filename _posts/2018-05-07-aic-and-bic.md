---
layout: post
title: AIC and BIC
date: 2018-05-07 23:01
author: baruuum
comments: true
categories: [AIC, BIC, Information Criteria, Information, Quant Stuff]
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos and errors._


Whenever several models are fitted to a dataset, the problem of model selection emerges. If the candidate models are nested the likelihood-ratio statistic or the F-test seems to be the preferred choice in the social science. For non-nested candidate models, on the other hand, the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC) are, by far, the most often used statistics. Both the AIC and the BIC are based on maximized likelihoods. Although it seems that the likelihood itself can be used to compare candidate models, it suffers from a similar problem as the R-squared statistics: namely, the the more parameters you fit to the data, the higher will be your likelihood, even if the model is wrong. Therefore, the likelihood itself is a bad choice for model comparison, except in the unique case where the ''dimensionality'' of the compared models are the same.

AIC and BIC have the form

$$  -2 \times (\text{maximized log-likelihood} - \text{penalty term}),  $$

so that they are often introduced as methods to penalize models for overfitting the data. Where these ''penalty terms'' come from is however not well-understood, except that the name ''information criterion'' suggests that they have something to do with information theory and that the BIC is related to Bayesian thinking. Some introduce AIC and BIC by using references to Ockham's razor. But as we will see below, they don't have anything to do with William of Ockham nor with the intrinsic value of parsimony in explanations.

So, what are AIC and BIC? As a student asked me this question, I realized that I had no good answer. So, I decided to look into the literature a little bit. This post of is the result of this search.

I start with defining the Kullback-Leibler Divergence between two distributions and thereafter derive the AIC. Then, I add a short discussion of the BIC.

In all what follows we assume that the data matrix $ x $ we observe (denoted in lower case for the sake of clarity in the following equations) is an iid sample from some ''true'' distribution $ p  $ which is not known. Also, we assume that all densities/distributions that are considered satisfy the regularity conditions, so that we can freely interchange the orders of limiting processes and can assume that there exists a unique value parameter value that maximizes the likelihood. The derivation of both AIC and BIC rely on truncations of Taylor series expansions. We will not derive the order of the error of any approximations that result these procedures and state them only when necessary, which are situations in which the order of the error is not well-known or when the error does not tend to zero with growing sample sizes.

<h3>Kullback-Leibler Divergence</h3>

Almost all so-called information criteria are derived from the <em>Kullback-Leibler (K-L) Divergence</em> of a model from the true distribution, and the AIC is no exception. Similarly to Shannon's entropy measure (see this [post](/blog/2018/entropy)), the K-L Divergence was derived from minimal assumptions regarding the information loss when we approximate one distribution with another. We have the following setup: suppose there exists a ''true'' distribution $ p  $ that is responsible for generating the data we are interested in. Of course, $ p  $ is not known to us. We set up a model of the data in terms of a <em>family</em> of distributions $ q(\cdot\vert\theta)=q_\theta  $, where $ \theta  $ is the vector of parameters of our model. For instance, we might let $ q_\theta  $ be the family of Normal distributions with mean $ \mu  $ and standard deviation $ \sigma  $, so that $ \theta = (\mu,\sigma)  $; different combinations of $ \mu  $ and $ \sigma  $ give rise to a different distributions, which all belong to the family of Normal distributions. Thus, our model $ q_\theta  $ specifies not any particular distribution, say $ \text{Normal}(1,.2)  $, but only $ \text{Normal}(\mu,\sigma)  $, where $ \theta  $ is left unspecified and will be estimated from the data. Notice also that, so far, we have no ''data'' at all. We have only the a distribution $ p  $ and a family of distributions $ q_\theta  $.

Next, fix $ \theta  $ at any particular value. Now $ q_\theta  $ <em>is</em> a particular distribution instead of a family of distributions. The K-L divergence <em>from</em> $ q_\theta  $ <em>to</em> $ p  $ is then given as

$$  D_{KL}(p\|q_\theta) = -\int p(x) \log \left(\frac{q_\theta(y)}{p(y)}\right)dy = \int p(y) \log \left(\frac{p(y)}{q_\theta(y)}\right)dy  $$

where $ y  $ is just the variable of integration and does not refer to any data (We could also have written $ \int_X p\log(p/q_\theta) d\mu  $, where $X$ is the set over which the distributions $p$ and $q_\theta$ are supported and where $\mu$ is a dominating measure.)
Notice that the K-L divergence might be rewritten as

$$ \begin{aligned}
D_{KL}(p\|q_\theta)&=E_p[\log(p)] - E_p[\log(q_\theta)] \\
&=-H(p) + H(p,q_\theta)
\end{aligned} $$

where we indicate by the subscript that the expectation is taken with respect to the ''true'' distribution $ p  $, and $ H(p)  $ and $ H(p,q_\theta)  $ are, respectively, the entropy of $ p  $ and the cross-entropy between $ p  $ and $ q_\theta  $. In information theoretic terms, $ D_{K}(p\\|q_\theta)  $ is the information loss that occurs when approximating $ p  $ (the truth) with $ q_\theta  $ (our model); heuristically, this might be understood as the ''directed'' distance from $ q_\theta  $ to $ p  $. It is clear that if $ p=q_\theta  $, $ D_{KL}(p\|q_\theta)=0  $. Also, by the Jensen's Inequality, we have that $ D_{KL}(p\|q_\theta)\ge 0  $. On the other hand, the K-L divergence is not symmetric, i.e., $ D_{KL}(p\|q_\theta)\ne D_{KL}(q_\theta\|p)  $, nor does it satisfy the triangle inequality, so that it is not a distance in a formal sense. As the divergence from $ p  $ to $ q_\theta  $ is not equal to that of $ q_\theta  $ to $ p  $, we have to be careful with the notation.

Now, suppose we don't know the truth $ p  $ but have a set of candidate models $ \mathcal Q=\\{q_1,q_2,...,q_m\\}  $ with which we want to approximate $ p  $, where we use the shorthand $ q_k  $ for $ q_k(\cdot \vert{\theta_k})  $. Let $ q_1,q_2\in \mathcal Q  $ be two candidate model and consider

$$  D_{KL}(p\|q_k) = E_p[\log p] - E_p[\log q_k],\quad k\in\{1,2\}.  $$

We find that the entropy of the ''truth,'' $ E_p[\log p] = \int p(y) \log p(y) dy  $ is a constant. So, when we are only interested in the <em>relative</em> K-L divergences between two candidate models from $ p  $, we get

$$  D_{KL}(p\|q_1)- D_{KL}(p\|q_2) = E_p[\log q_2] - E_p[\log q_1].  $$

In words, ''truth'' is subtracted out and only the cross-entropy remains. Hence, we do not need to know the entropy of the true distribution in order to calculate the relative divergence from the truth of the candidate models: the quantity $ E_p[\log q(x\vert\theta_k)]  $ is sufficient to compare them.

Still, we need to know $ p  $ in order to calculate $ E_p[\log q(\cdot\vert\theta_k)]  $ as it is an average over the distribution $ p  $. But, as assumed, $ p  $ is not known. Hence, we start with a different question: namely, what is the value $ \theta_0\in \Theta  $, $ \Theta  $ denoting the parameter space, such that the K-L divergence of $ q(\cdot \vert \theta)  $ to $ p  $ is minimized? That is, we are seeking the value $ \theta_0  $ that solves the problem

$$  \min_{\theta\in \Theta}D_{KL}(p\|q_{\theta})=\min_{\theta\in \Theta}\int p(y) \log\left(\frac{p(y)}{q(y\vert\theta)}\right) dy.  $$

If there exists such a $ \theta_0  $, it has to satisfy the first-order condition

$$  \frac{\partial}{\partial \theta} \int p(y)\log q(y\vert\theta)dy = E_p\left[\frac{\partial}{\partial \theta}\log q(y\vert\theta_0)\right]=0.  $$

This equation looks familiar...Indeed, if we have an iid sample of size $ n  $ from $ p  $ (not $ q_\theta  $) and we are modeling it with $ q_\theta  $, then the maximum likelihood estimate of $ \theta  $, which we denote by $ \tilde\theta  $ must satisfy the first-order condition

$$  \frac{1}{n}\left[\sum_{i=1}^n \frac{\partial}{\partial \theta}\log q(x_i\vert\tilde\theta)\right] = 0.  $$

Further, by the Strong Law of Large Numbers, it follows that

$$  \frac{1}{n}\left[\sum_{i=1}^n \frac{\partial}{\partial \theta}\log q(x_i\vert\tilde\theta)\right] \longrightarrow E_p[(\partial/\partial \theta)q(y\vert\tilde\theta)] \qquad\text{almost surely}. $$

In other words, the value to which the maximum likelihood estimator $ \tilde\theta  $ converges as the sample becomes large is exactly the value of $ \theta  $ that solves the minimum K-L divergence problem. As this value is unique by assumption (this follows from the regularity conditions), it has to be the case that $ \tilde\theta \rightarrow \theta_0  $ as $ n \rightarrow \infty  $. That is, among the family of distributions $ q(\cdot\vert\theta)  $, $ \theta\in \Theta  $, the distribution specified by the maximum likelihood estimator $ \tilde\theta  $, $ q(\cdot\vert\tilde\theta)  $, minimizes the K-L divergence of $ q_\theta  $ from $ p  $.

Now, suppose that $ p  $ is indeed within the family of distributions $ q(\cdot\vert\theta)  $. That is, there exists a value $ \theta_0\in \Theta  $ such that $ p=q(\cdot\vert\theta_0)=q_{\theta_0}  $. In such a situation we would have $ D_{KL}(p\|q_{\theta_0}) = 0  $. However, we have an additional problem: namely, we have to <em>estimate</em> $ \theta_0  $ with our data $ x  $ using $ \tilde\theta(x)  $, where the parentheses behind $ \tilde\theta  $ are added just to emphasize that $ \tilde\theta  $ depends on $ x  $. For each particular sample $ x  $, we will obtain different estimates $ \tilde\theta(x)  $ and, accordingly, different estimates of the cross entropy. Actually, the probability that $ D_{KL}(p\|q_{\tilde\theta}) = 0  $ for any specific sample is zero. Thus, we have to base our model selection procedure on a different quantity.

A reasonable choice would be to choose the model that is closest to ''truth'' <em>in expectation</em>. That is, we want the model that minimizes

$$  Q=E_{x\sim p}\Big[D_{KL}\big(p\|q(\cdot \vert\tilde \theta(x)\big)\Big]  $$

where the subscript $x\sim p$ means that the expectation is taken over repeated samples of $ x  $ from $ p  $. We note that

$$ \begin{aligned}
E_{x\sim p}\Big[D_{KL}\big(p\|q(\cdot \vert\tilde \theta(x)\big)\Big] &= E_{x\sim p}\left[\int p(y) \log p(y) dy\right] - E_{x\sim p}\left[ \int p(y)\log q\big(y\big\vert\tilde\theta(x)\big)dy \right] \\
&=\text{ Constant } - E_{x\sim p}\left[E_{y\sim p}\Big[\log q\big(y \big\vert \tilde\theta(x)\big)\Big]\right],
\end{aligned} $$

where both expectations, as indicated by the subscript, are taken with respect to truth $ p  $. So, minimizing the expected K-L divergence will be equivalent to maximizing the double expectation:

$$  Q = E_{x\sim p}\left[E_{y\sim p}\Big[\log q\big(y \big\vert \tilde\theta(x)\big)\Big]\right]= \int \left[\int \log q(y\vert\tilde\theta(x)) p(y) dy\right]p(x) dx.  $$

In intuitive terms, $ Q  $ can be understood as a quantity that we obtain when we repeat the following process an infinite number of times and take the average of the resulting quantities: sample a dataset $x$ from $p$ and maximize $ \log q  $ to obtain $ \tilde\theta(x)$; thereafter sample another dataset $y$ from $p$ and evaluate $q_\theta(y\,\vert\,\tilde\theta(x))$ on $y$. It follows that $ Q  $ can be interpreted as a kind of ''expected out-of-sample log-likelihood.''

$$  Q = E_p\left[E_p\left[\log q(y\vert\tilde\theta(x))\Big\vert x\right]\right]= E_p\left[E_p\left[\log\mathcal L(\tilde\theta(x))\Big\vert x\right]\right]= E_p\left[\log \mathcal L(\tilde\theta)\right].  $$

AIC is, in essence, an estimator of $ Q  $ multiplied by $ -2  $ for historical reasons (as most goodness-of-fit statistics are on the deviance scale). The interpretation as an expected out-of-sample log-likelihood is the connection that links the AIC to cross-validation.

So far, we made clear that AIC has nothing to do with Ockham's Razor, but is based on the K-L divergence between truth $ p  $ and our model $ q_\theta  $. The maximized log-likelihood from the data $ x  $, $ \mathcal L(\tilde\theta(x))  $, might seem to be a good approximation to $ Q  $, the expected out-of-sample log-likelihood. Unfortunately, this is not the case in that $ \mathcal L(\tilde\theta)  $ will be upwardly biased. The seminal insight of Akaike in deriving the AIC was to show that

$$  \mathcal L(\tilde\theta(x)) - E_p\left[\log \mathcal L(\tilde\theta)\right] \approx \text{dim}(\theta),  $$

in ''good models, '' where $ \text{dim}(\theta)  $ is the dimensionality of $ \theta  $, which can be roughly interpreted as the number of parameters in $ \theta  $. What we mean by ''good models'' will become clear shortly.

<h3>Derivarion of AIC</h3>

To derive the AIC, we first Taylor expand $ h(\tilde\theta) = \log q(y\vert\tilde\theta(x))  $ around $ \theta_0 $, the value in $\Theta$ that maximizes the likelihood. For clarity, we omit the dependence of $ \tilde\theta  $ on the data $ x  $ in the notation.

$$  h(\tilde\theta) \approx h(\theta_0) + h'(\theta_0)^\top (\tilde\theta - \theta_0) + \frac{1}{2}(\tilde\theta - \theta_0)^{\top}h''(\theta_0) (\tilde\theta - \theta_0),  $$

where

$$ \begin{aligned}h'(\theta_0) &= \left[\frac{\partial}{\partial \theta} \log q(y \vert\theta) \right]_{\theta=\theta_0} \\
h''(\theta_0)&= \left[\frac{\partial^2}{\partial \theta \partial \theta^{\top}} \log q(y \vert\theta) \right]_{\theta=\theta_0}
\end{aligned} $$

are, respectively, the gradient vector and Hessian matrix evaluated at $ \theta_0  $. Taking expectations with respect to $ y  $ (which has distribution $ p  $ and is independent of $ x  $), we obtain

$$ \begin{aligned}E_y[h(\tilde\theta)] &\approx E_{y\sim p}[h(\theta_0)] + E_{y\sim p}\left[\frac{1}{2}(\tilde\theta - \theta_0)^{\top}h''(\theta_0) (\tilde\theta - \theta_0)\right]\\
&=E_{y\sim p}[\log q(y\vert\theta_0)] - \frac{1}{2}(\tilde\theta - \theta_0)^{\top}I_p(\theta_0) (\tilde\theta - \theta_0)
\end{aligned} $$

since the linear term vanishes as $ E_{y\sim p}[h'(\theta_0)] = E_p[(\partial /\partial \theta)\log q(y\vert\theta_0)] = 0  $ by the first order condition and where

$$  I(\theta_0) = E_p[-h''(\theta_0)] = E_p\left[-\frac{\partial^2}{\partial \theta\partial \theta^\top}\log q(y\vert\theta)\right]_{\theta=\theta_0},  $$

i.e., the expected negative Hessian evaluated at $ \theta_0  $. Again, the expectation is with respect to $ p  $ not $ q  $, and so $ I(\theta_0)  $ matrix is, in general, not the Fisher information matrix ( $ I(\theta_0)  $ is equal to the Fisher information only if $ p=q  $).

As the next step, we take expectations with respect to $ x  $. By noting that only $ \tilde\theta  $ in the expansion depends on $ x  $, that $ x  $ and $ y  $ are independent, and that

$$ \begin{aligned}
(\tilde\theta - \theta_0)^{\top}I(\theta_0) (\tilde\theta - \theta_0) = \text{tr}\left[(\tilde\theta - \theta_0)^{\top}I(\theta_0) (\tilde\theta - \theta_0)\right] \\
=\text{tr}\left[I(\theta_0) (\tilde\theta - \theta_0)(\tilde\theta - \theta_0)^{\top}\right]
\end{aligned} $$

by the cyclic property of matrix trace, we obtain

$$  E_{x\sim p}\left[E_{y\sim p}\left[q\big(y\big\vert\tilde\theta(x)\big)\right]\right] \approx E_{y\sim p}[\log q(y\vert\theta_0)] - \frac{1}{2}\text{tr}\left\{I(\theta_0)E_{x\sim p}\left[(\tilde\theta - \theta_0)(\tilde\theta - \theta_0)^{\top}\right]\right\}.  $$

The left-hand side is $ Q  $ and $ E_{x\sim p}\{[\tilde\theta(x)-\theta_0][\tilde\theta(x)-\theta_0]^\top\}=\Sigma,  $ the asymptotic covariance matrix of the MLE. Thus, we have the result

$$  Q \approx E_{y\sim p}\left[\log q(y\vert\theta_0)\right] -\frac{1}{2}\text{tr}\left[ I(\theta_0)\Sigma\right].  $$

The problem that remains is to estimate the terms on the right hand-side of the equation from the data. To derive an estimator, we Taylor-expand $ \log q(x\vert\theta_0)  $, this time around $ \tilde\theta(x)  $, and take expectations with regard to $ x  $ of both sides. This results in

$$  E_{x\sim p}\left[\log q(x\vert\theta_0)\right] \approx E_{x\sim p}\left[\log q(x\vert\tilde\theta)\right] - \frac{1}{2}\text{tr}\left\{E_{x\sim p}\left[ I(\tilde\theta)(\theta_0-\tilde\theta)(\theta_0-\tilde\theta)^\top\right]\right\}  $$

where $ I(\tilde\theta)  $ is the negative Hessian, this time evaluated at the MLE. The derivation is almost identical to the first expansion from above. Next, we approximate $ I(\tilde\theta)  $ with $ I(\theta_0)  $, which does not depend on $ x  $, to get

$$ \begin{aligned}
E_{x\sim p}\left[\log q(x\vert\theta_0)\right] &\approx E_{x\sim p}\left[\log q(x\vert\tilde\theta)\right] - \frac{1}{2}\text{tr}\left\{ I(\tilde\theta)E_{x\sim p}\left[(\theta_0-\tilde\theta)(\theta_0-\tilde\theta)^\top\right]\right\}\\
&=E_{x\sim p}\left[\log q(x\vert\tilde\theta)\right] - \frac{1}{2}\text{tr}\left\{I(\theta_0)\Sigma\right\}.
\end{aligned} $$

By observing that $ E_{x\sim p}[\log q(x\vert\theta_0)] = E_{y\sim p}[\log q(y\vert\theta_0)]  $, as both $ x  $ and $ y  $ come from the same distribution $ p  $, and substituting this result into the expression for $ Q  $, we obtain

$$  Q \approx E_{x\sim p}[\log q(x\vert\tilde\theta(x))] - \text{tr}\left\{I(\theta_0)\Sigma\right\}.  $$

Further, using the large sample result

$$  I(\theta_0)\Sigma = J(\theta_0)I(\theta_0)^{-1},  $$

where

$$  J(\theta_0)= E_p\left[\left[\frac{\partial}{\partial \theta} \log q(x\vert\theta)\right]_{\theta=\theta_0}\left[\frac{\partial}{\partial \theta} \log q(x\vert\theta)\right]_{\theta=\theta_0}^{\top}\right],  $$

we rewrite the $ Q  $ as

$$  Q\approx E_{x\sim p}\left[\log q(x\vert\tilde\theta(x))\right] - \text{tr}\{J(\theta_0)I(\theta)^{-1}\}.  $$

(Small note: the result $ I(\theta_0)\Sigma = J(\theta_0)I(\theta_0)^{-1}  $ is what leads to ''robust'' variance estimators in the MLE context. The result can be easily derived by using Taylor-expansion of the likelihood equations around $ \theta_0  $ and taking expectations with respect to $ p  $.) Since the maximized log-likelihood is an unbiased estimator for itself, a natural estimator for $ Q  $ is thus

$$  \hat Q_{TIC} = \log q(x\vert\tilde\theta(x)) - \widehat{\text{tr}}\{J(\theta_0)I(\theta)^{-1}\}.  $$

which is the definition of the TIC (Takeuchi's Information Criterion).

To derive the AIC, we further make the assumption that our model $ q  $ is a ''good model,'' i.e., that it is close to $ p  $. If so, we have $ J(\theta_0)\approx I(\theta_0) \approx \mathcal I(\theta_0)  $, by the Fisher Information Equality (see discussion on calculating standard errors of [this post](/blog/2017/EM/)), where $ \mathcal I(\theta_0)  $ is the Fisher Information. Thus for ``good models,'' $ J(\theta_0)I(\theta_0)^{-1} \approx \text{I}_d  $, the $ d  $-dimensional identity matrix. It follows that that $ \text{tr}(\text{I}_d)=d  $ and the quantity

$$  \hat Q = \log(q\vert\tilde\theta(x)) - d  $$

gives us a simple approximation to $ E_{x\sim p}\left[E_{y\sim p}\Big[\log q\big(y \big\vert \tilde\theta(x)\big)\Big]\right]  $, the maximization of which minimizes the expected K-L divergence from the $p$. Lastly, out of ''historical'' reasons, the AIC is defined by multiplying $ \hat Q  $ by $ -2  $, which yields, the AIC:

$$  AIC = -2 \log q(x\vert\tilde\theta(x)) + 2d.  $$

To summarize, the AIC is an approximation to the expected cross entropy between the true distribution $ p  $ and a set of models $ \mathcal Q  $. It turns out that the difference in the expected cross entropy between two models and the true distribution is equal to the <em>relative</em> Kullback-Leibler divergence between the models and $ p  $. Using several Taylor expansions, we derive at an estimator for the expected cross entropy between $ p  $ and our model. This estimator is given by the log-likelihood evaluated at the MLE, $ \log q(x\vert\tilde\theta(x))  $ minus a ''correction term,'' which is equal to $ \text{tr}(J(\theta_0)I(\theta_0)^{-1})  $, where $ J(\theta_0)  $ is the expected outer product of the score vector and where $ I(\theta_0)^{-1}  $ is the expected negative Hessian both evaluated at $ \theta_0  $ and where the expectation is over the true distribution $ p  $. For ''good models'' we have $ J(\theta_0) \approx I(\theta_0)  $, so that we can estimate the correction term by the dimensionality of $ \theta  $, $ \text{dim}(\theta) = d  $, i.e., the number of parameters in the model. Multiplying $ \hat Q = \log q(x\vert\tilde\theta(x)) - d  $ by $ -2  $ gives rise to the AIC. As larger values of $ \hat Q  $ indicate a smaller divergence from the ''truth'' $ p  $ in the K-L sense, models with smaller values of AIC are preferred.

<h3>Derivation of the BIC</h3>

In some sense, the Bayesian Information Criterion (BIC) is a double-misnomer: it is not a information criterion in the sense of not being based on Kullback-Leibler divergence, and it is not truly Bayesian, as it uses the maximum likelihood estimate of the parameter, rather than averaging over the posterior distribution. Anyways, the basic setup is as follows: We, again, have a set of candidate models $ \mathcal Q=\\{q_1,q_2,...,q_m\\}  $ and we want a procedure to selects the best model. As before, choosing the model with the highest likelihood will not work as ''larger'' models will have higher likelihoods, even if they are the ''wrong models.'' Instead, the BIC tries to solve the model selection problem by approximating the logarithm of the so-called _marginal likelihood_ or _model evidence_.

Let $ \\{\theta_1,...,\theta_m\\}  $ be a parameter vectors specified by the models in $ \mathcal Q  $. The set $ \mathcal Q  $ might be conceived as competing hypotheses, each of them encoded into a statistical model. We denote the observed data matrix by $ x  $, and are interested in the probability of a model $ q_k  $ being the correct model given the data. By Bayes' Theorem, this probability is

$$  p(q_k\vert x) =\frac{p(x\vert q_k)p(q_k)}{p(x)},  $$

where $ p(q_k)  $ is the prior probability we put on $ q_k  $ being the right model and $ p(x\vert q_k)  $ is the probability of the data under the model $ q_k  $. So, in this set up we are treating $ p  $ as a generic term of a probability distribution and **not** the true model. The ''model evidence'' or the ''marginal likelihood'' is the likelihood of the data after the parameters of the model are ''marginalized out.'' Formally,

$$  p(x\vert q_k) = \int p(x, \theta_k \vert q_k) d\theta_k = \int p(x\vert\theta_k, q_k)p(\theta_k\vert q_k)d\theta_k  $$

where the $ p(x\vert\theta_k,q_k)  $ is the likelihood of the data and $ p(\theta_k\vert q_k)  $ is the prior distribution of the parameter vector under model $ q_k  $.

The marginal likelihood plays an important role in model selection and the ratio of two marginal likelihoods is known as the _Bayes factor_. For instance, when we have two competing models $ q_1  $ and $ q_2  $, the Bayes factor of $ q_1  $ over $ q_2  $ is defined as

$$  B(q_1,q_2) = \frac{p(x\vert q_1)}{p(x\vert q_2)}.  $$

Note that multiplying $ B(q_1,q_2)  $ by the ratio of the prior odds $ p(q_1)/p(q_2)  $ gives the posterior odds:

$$  \frac{p(q_1\vert x)}{p(q_2\vert x)} = B(q_1,q_2)\frac{p(q_1)}{p(q_2)}  $$

as the normalizing constant for both the numerator and the denominator are the same.

The BIC can be understood as an asymptotic approximation to the log of the marginal likelihood with a specific prior for the model parameters. For ease of exposition, we consider a single generic model $ q\in \mathcal Q  $. Also, let $ \theta  $ be the parameter vector associated with $ q  $. Keeping in mind that $ p(x\vert\theta)=p(x\vert\theta, q)  $ and $ p(\theta)=p(\theta\vert q)  $ we drop to reference to $ q  $ for clarity in what follows.

Now, define

$$  g(\theta) =\log \Big[p(x\vert\theta) p(\theta)\Big].  $$

Note that $ p(x\vert\theta)p(\theta)\propto p(\theta\vert x)  $, the posterior distribution of $ \theta$. The marginal likelihood is the normalization constant of the posterior and, thus, can be written as $ p(x) = \int e^{g(\theta)} d\theta  $. Let $ \hat\theta= \text{argmax}_\theta p(x\vert\theta)p(\theta)  $, the value of $ \theta  $ that maximizes the posterior of $ \theta  $. $ \hat\theta  $ is simply the posterior mode and is also called the <em>maximum a posteriori</em> (MAP) estimate of $ \theta  $.

To approximate the marginal likelihood, $\int \exp[g(\theta)]d\theta$, we use the <em>Laplace approximation</em>. Consider the Taylor expansion of $ g  $ near the point $ \hat\theta  $,

$$  g(\theta) \approx g(\hat\theta) + g'(\hat\theta)^\top (\theta-\hat\theta) + \frac{1}{2} (\theta - \hat\theta)^\top g''(\hat\theta)(\theta-\hat\theta)  $$

where $ g'(\hat\theta)  $ is the gradient vector and $ g''(\hat\theta)  $ is the Hessian matrix evaluated at $ \hat\theta  $. As $ \theta  $ achieves its maximum at $ \hat\theta  $, the linear term has to be zero. Thus, we have

$$  g(\theta) \approx g(\hat\theta) + \frac{1}{2} (\theta - \hat\theta)^\top g''(\hat\theta)(\theta-\hat\theta).  $$

By exponentiating both sides and integrating, we can approximate the marginal likelihood by

$$  p(x) \approx \exp\big[g(\hat\theta)\big] \int \exp\Big[\frac{1}{2}(\theta-\hat\theta)^\top g''(\hat\theta)(\theta-\hat\theta)\Big]d\theta.  $$

The integrand is the kernel of a multivariate Normal distribution with mean vector $ \hat\theta  $ and covariance matrix $ -[g''(\hat\theta)]^{-1}  $. So, by letting $ I(\hat\theta) = -g''(\hat\theta)  $, we have

$$  \int\exp\Big[\frac{1}{2}(\theta-\hat\theta)^\top g''(\hat\theta)(\theta-\hat\theta)\Big]d\theta = (2\pi)^{d/2}(\det I(\hat\theta))^{-1/2}  $$

as $ \det A^{-1} = [\det A]^{-1}  $, and where $ d  $ is the number of parameters in the model. Plugging these results back into our approximation of $ p(x)  $ and taking logs on both sides, we obtain

$$  \log p(x) \approx \log p(x\vert\hat\theta) + \log p(\hat\theta) + \frac{d}{2}\log(2\pi) - \frac{1}{2}\log\det I(\hat\theta).  $$

In large samples, the posterior becomes dominated by the likelihood (recall that the prior is weighted equal to ''one observation,'' while the likelihood consists of $ n  $).This implies that the MAP estimator, $ \hat\theta  $, becomes arbitrarily close to the MLE, $ \tilde\theta  $, and that $ I(\hat\theta)  $ approaches the Fisher Information matrix $ \mathcal I(\tilde\theta)  $. Thus, we approximate $ \hat\theta  $ by $ \tilde\theta  $ and $ I(\hat\theta)  $ by $ \mathcal I(\tilde\theta)  $. In addition, noting that

$$  (1/2)\log[\det I(\hat\theta)] \approx (1/2)\log [\det(n \mathcal I_i(\tilde\theta))] = (1/2) \log[n^d \det(\mathcal I_i(\tilde\theta))],  $$

where $ \mathcal I_i(\tilde\theta)  $ is the Fisher information for a single observation, the last equation can be rewritten as

$$  \log p(x) \approx \log p(x\vert\tilde\theta) + \log p(\tilde\theta) + \frac{d}{2}\log(2\pi) - \frac{d}{2}\log n - \frac{1}{2}\log\det \mathcal I_i(\tilde\theta)  $$

Notice that $ \log p(\tilde\theta)  $, $ (d/2)\log(2\pi)  $, and $ \mathcal I_i(\tilde\theta)  $ are constants, not depending on $ n  $. So,

$$  \log p(x)= \log p(x\vert\tilde\theta) - \frac{d}{2}\log n + O(1).  $$

Our estimator to approximate the log marginal likelihood is therefore

$$  \hat Q = \log p(x\vert\tilde\theta) - \frac{d}{2}\log n,  $$

which is, similarly to the AIC, the log-likelihood evaluated at the MLE minus a ''correction term.'' The error $ O(1)  $ implies that the bias does not go away even in infinite samples. Yet, with an appropriate choice of $ p(\theta)  $, this problem might be overcome. For example, if the support of $ \theta  $ is $ \mathbb R^{d}  $, then by choosing $ p(\theta) \sim \text{Multivariate Normal}_d(\tilde\theta, \mathcal I_i(\tilde\theta)^{-1})  $, we obtain

$$  \log p(\tilde\theta) = -\frac{d}{2}\log(2\pi) + \frac{1}{2}\log \det\mathcal I_i(\tilde\theta).  $$

Thus, for this choice of prior, we have

$$  \log p(x) =\log p(x\vert\tilde\theta) - \frac{d}{2}\log n + O(n^{-1/2}),  $$

where $ O(n^{-1/2})  $ is the error rate introduced by approximating and $ \hat\theta  $ by $ \tilde\theta  $ and $ I(\hat\theta)  $ by $ \mathcal I(\tilde\theta  $) (the error rate of the Laplace approximation is of order $ O(n^{-1})  $ and is thus less than the second approximation that we have used). Hence, with an appropriately chosen prior, $ \hat Q  $ converges to the log-marginal-likelihood as $ n  $ tends to infinity. Again, for ''historical reasons,'' BIC is defined as $ -2  $ times the approximated log marginal likelihood, i.e.,

$$  BIC = -2\log p(x\vert\tilde\theta) + d\log n.  $$

As higher values of $ \hat Q  $ indicate higher marginal likelihoods, lower values of BIC are indications of higher posterior odds in favor of the model.

To summarize, the BIC is an approximation to the log-marginal-likelihood using a specific prior on the parameter vector (and multiplied by $-2$ for historical reasons). It's asymptotic consistency depends on the prior of the parameter vector; but when using the right one, e.g., Normal distribution when the support of $ \theta  $ is $ \mathbb R^d  $, BIC will be a consistent estimator of $ \log p(x\vert q)  $. Similarly to the AIC, lower values of BIC indicate "better" models. From the definition of AIC and BIC, it is clear that the BIC will prefer smaller models as long as $ n\ge 8  $. When the number of observations is large, conclusions based on AIC and BIC might be markedly different.

Another point that becomes clear from the derivation of AIC and BIC is that these measures cannot be compared across different data. In the case of AIC, the "truth" that has generated $ x_1  $ might not be the same as that which has generated $ x_2  $ and subtracting two AIC values will not give the relative divergence from the _same_ distribution. A similar reasoning applies also to the BIC. A little bit more ambiguous situation arises when one model is fitted on $ x  $ and another one on $ \log x  $ (for example, a regression where $ x  $ is the outcome). In principle, it should be possible to compare AIC values across these models, when the values are appropriately ''scaled''. Even more complicated is a situation in which we have hierarchical models. This is because the dimensionality of the parameter vector $ \theta  $ is not clear. For example, in a random effects model with two regressors (including the constant term) and a random intercept that varies over $ J  $ individuals, is the number of parameters $ 4  $ (the constant, the regression slope, the residual variance, and the variance of the random effect) or is it $ 3+(J-1)  $ (constant, regression regression slope, residual variance, and the $ J-1  $ random effects), or is it somewhere between $ 4  $ and $ 3+(J-1)  $? This problem leads to the concept of ''effective number of parameters'' of a model that is used in the derivation of the DIC (Deviance Information Criterion) and the WAIC (Watanabe-Akaike Information Criterion).
