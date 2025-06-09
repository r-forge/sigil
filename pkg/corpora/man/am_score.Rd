\name{am.score}
\alias{am.score}
\alias{builtin.am}
\title{Compute association scores for collocation analysis (corpora)}

\description{

  Compute a wide range of established association measures (Evert 2004, 2008) 
  for collocation analysis based on the co-occurrence frequency of two lexical
  items as well as their individual \dQuote{marginal} frequencies.
  Frequency data can be provided in a number of different ways for convenience,
  and user-defined measures can be specified as mathematical formulae using the
  flexible notation introduced by Evert (2004).
  
}

\usage{

am.score(w1, w2, f, f1, f2, N, measure,
         span.size=1, log=FALSE, labels=FALSE,
         conf.level=.95, p.adjust=TRUE,
         param=list())

builtin.am

}

\arguments{

  \item{w1}{a character vector specifying the first lexical item of each co-occurrence pair. If of length 1, it is assumed to be the same for all pairs. May be omitted if \code{f1} is an unnamed vector parallel to \code{f}.}
  \item{w2}{a character vector specifying the second lexical item of each co-occurrence pair. If of length 1, it is assumed to be the same for all pairs. May be omitted if \code{f2} is an unnamed vector parallel to \code{f}.}
  \item{f}{a numeric vector specifying the co-occurrence frequency of each pair of lexical items \eqn{(w_1, w_2)}}
  \item{f1}{a numeric vector of first marginal frequencies \eqn{f_1}, labelled with the unique strings occurring in \code{w1}. Alternatively, an unlabelled numeric vector of the same length as \code{f} or of length 1 (implying that all pairs have the same first marginal \eqn{f_1}).}
  \item{f2}{a numeric vector of second marginal frequencies \eqn{f_2}, labelled with the unique strings occurring in \code{w2}. Alternatively, an unlabelled numeric vector of the same length as \code{f} or of length 1 (implying that all pairs have the same first marginal \eqn{f_2}).}
  \item{N}{a numeric vector specifying the sample size \eqn{N} underlying the contingency tables. It is normally a single value, which is the same for all co-occurrence pair, but a vector of the same length as \code{f} may also be specified.}
  \item{measure}{either a string specifying one of the built-in association measures or a function implementing a user-defined measure (see \dQuote{Details} below)}
  \item{span.size}{for surface co-occurrence, specify the total size of each collocational span in tokens (e.g. \code{span.size=6} for a L3/R3 span) in order to adjust first marginal frequencies as explained by Evert (2008: Sec. 4.1 + 5.1). No adjustment should be specified if contingency tables have already been obtained according to Evert (2004: 68-70).}
  \item{log}{if TRUE, apply a signed log transformation \eqn{h(x) = \mathop{\text{sgn}}(x)\cdot \log_2(|x| + 1)}{h(x) = sgn(x) * log2(|x| + 1)} to the computed association scores in order to deskew their distribution.}
  \item{labels}{if TRUE, add labels to the returned vector of assocation scores, constructed from \code{w1} and \code{w2}}
  \item{conf.level}{the desired confidence level for association measures based on confidence intervals such as \code{LRC} (defaults to 95\%)}
  \item{p.adjust}{if \code{TRUE}, apply a Bonferroni correction in order to control the family-wise error rate across all statistical tests carried out across the data set. This parameter is used by measures based on confidence intervals and measures that return (log) p-values.}
  \item{param}{a list of additional parameters passed to the `param` argument of the selected association measure, with `conf.level` and `p.adjust` added automatically.}
}

\value{

  A numeric vector of the same length as \code{f} containing the desired association scores for all co-occurrence pairs. The vector is labelled with the corresponding co-occurrence pairs if \code{label=TRUE}.
  
  For all standard association measures, larger values indicate higher collocability. There is usually no meaningful scale  for interpretation of the scores, though many measures return positive scores for higher than expected co-occurrence frequency and negative scores for lower than expected co-occurrence frequency. 

  P-values are returned on a negative log-10 scale (\eqn{-\log_{10} p}), so a score of 2 corresponds to \eqn{p = .01} and a score of 3 to \eqn{p = .001}. The score for \eqn{p = .05} is approximately 1.3. These cutoffs can be used to implement a significance filter for collocation candidates.
  
}

\details{

  This function computes a range of standard association measures for collocation analysis, as well as user-defined measures specified by mathematical formulae. It assumes familiarity with the basic approaches of collocation analysis as explained in Evert (2008) or in Unit 4 of the SIGIL course at \url{https://SIGIL.R-Forge.R-Project.org/}, as well as with the notation for contingency tables introduced there.
  
  While the main purpose of the function is to determine the most strongly associated items in a large data set of co-occurrence pairs, it can also be applied to the collocates of a single node word. In this case, \code{w1} and \code{f1} only consist of a single item (the node word and its marginal frequency). Keep in mind that marginal frequencies \code{f2} still need to be obtained from the entire corpus, not just from co-occurrences with the node.
  
  \subsection{Frequency data}{
  
    Frequency signatures for all co-occurrence pairs \eqn{(w_1, w_2)} must be provided by the user in terms of a vector of co-occurrence frequencies \code{f}, as well as vectors of marginal frequencies \code{f1} (for the lexical items \eqn{w_1}) and \code{f2} (for the lexical items \eqn{w_2}). The underlying sample size \code{N} of the contingency tables also has to be specified, which is usually the same for all co-occurrence pairs.
  
    Marginal frequencies are normally given as numeric vectors labelled with the distinct strings from \code{w1} (for \code{f1}) and the distinct strings from \code{w2} (for \code{f2}), respectively. These vectors are used to look up the appropriate marginal frequencies for each co-occurrence item based on \code{w1} and \code{w2}. Alternatively, \code{f1} and \code{f2} can be unlabelled vectors of the same length as \code{f} (or of length 1), in which case \code{w1} and \code{w2} may be omitted.

    Because of these different invocations, it is recommeded to always use named arguments to \code{am.score()}.
    
    Contingency tables are automatically derived from the frequency signatures internally, and can be accessed by user-defined measures using the flexible and convenient notation of Evert (2004, 2008).

  }

  \subsection{Types of co-occurrence}{
  
    Evert (2008: Sec. 3) distinguishes three types of co-occurrence, which have implications on how frequency signatures and contingency tables have to be constructed. While reading the full description in Evert (2008) is strongly encouraged, concise summaries are provided here for reference.
    
    \bold{Syntactic co-occurrence} only considers pairs of lexical items in a specific syntactic relation to each other, such as adjective-noun or verb-object. Frequency data are obtained by extracting all instances of this syntactic relation from a corpus, resulting in a list of pair tokens. The sample size \eqn{N} is the total number of such pair tokens, while co-occurrency frequency \eqn{f} is obtained by counting the pairs, and marginal frequencies \eqn{f_1, f_2} by counting lexical items in the first or second element of the pairs, respectively. For a \dQuote{word sketch} (or \dQuote{word profile}) collecting multiple syntactic relations, the collocation analysis for each relation has to be carried out separately.
  
    \bold{Textual co-occurrence} refers to co-occurrence within text segments such as sentences, paragraphs, tweets, or chapters. In this case, the sample size \eqn{N} corresponds to the total number of relevant text segments, and \eqn{f_1, f_2, f} are \dQuote{document frequencies}, i.e. the number of text segments in which \eqn{w_1, w_2}, or both lexical items occur.
    
    \bold{Surface co-occurrence} usually considers \eqn{w_1} as the \dQuote{node} of the analysis and looks for collocates \eqn{w_2} that occur within a certain number of tokens around the node (the \dQuote{collocational span}). The span can be symmetric, e.g. L3/R3 (three tokens to the left and right of the node) or asymmetric, e.g. L0/R5 (five tokens to the right of the node only). The span size \eqn{k} is the total number of tokens in each span (6 for L3/R3, and 5 for L0/R5). Evert (2008) recommends to obtain frequency signatures as follows: \eqn{f} from the co-occurrence counts within collocational spans, \eqn{f_1} as the overall corpus frequency of each distinct node, \eqn{f_2} as the overall corpus frequency of each distinct collocate, and sample size \eqn{N} as the total number of tokens in the corpus. Contingency tables and expected frequencies then need to be adjusted for the span size (Evert 2008: Sec. 4.1 + 5.1), which can easily be achieved with the argument \code{span.size=k}.
  
  }

  \subsection{Built-in association measures}{
    
    The \code{am.score()} function includes a wide selection of built-in association measures, which can be selected via their name in the \code{measure} argument. Please refer to Evert (2004, 2008) or \url{http://www.collocations.de/AM/} for mathematical details on these measures and their properties as well as full equations. Formulae shown below use the flexible notation for contingency tables introduced there.
    
    Some measures take additional paramters specified in the \code{param} argument. Measures that return (log) p-values apply a Bonferroni correction if \code{p.adjust=TRUE} (the default). Measures based on confidence intervals use the confidence level specified in \code{conf.level}, which is also adjusted by the Bonferrroni correction.
  
      \describe{
        \item{\code{MI}}{
        (Pointwise) \bold{mutual information}, the binary logarithm of the ratio between observed and expected co-occurrence frequency: \deqn{
          \log_2 \frac{O_{11}}{E_{11}}
        }{ log2(O11 / E11) }
        Pointwise MI has a very strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). It should usually be combined with frequency thresholds on \eqn{f_1} and \eqn{f_2}.
      }
    
      \item{\code{MI.k}}{
        A heuristic variant of pointwise \bold{mutual information} intended to counteract the low-frequency bias by raising \eqn{O_{11}}{O11} to power \eqn{k}: \deqn{
          \log_2 \frac{(O_{11})^k}{E_{11}}
        }{ log2(O11^k / E11) }
        The exponent \eqn{k} can be specified as a user parameter (e.g. \code{param=list(k=3)} for the popular MI3 measure) and defaults to \eqn{k=2} otherwise.
      }
      
      \item{\code{G2}}{
        The \eqn{G^2} statistic of a likelihood ratio test for independence of rows and columns in a contingency table, which is very popular in computational linguistics under the name \bold{log-likelihood}: \deqn{
          \pm 2 \left( \sum_{ij} O_{ij}\cdot \log \frac{O_{ij}}{E_{ij}} \right)
        }{ ± 2 * ( SUM[ij] Oij * log(Oij / Eij) ) }
        This implementation computes \emph{signed} association scores, which are negative iff \eqn{O_{11} < E_{11}}{O11 < E11}.
        Log-likelihood has a strong bias towards high co-occurrence frequency and often produces a highly skewed distribution of scores. It may therefore be advisable to combine it with an additional log transformation (\code{log=TRUE}).
      }
      
      \item{\code{G2.pv}}{
        The \bold{p-values} corresponding to the \eqn{G^2} scores of the likelihood ratio test. In order to achieve sensible scaling and ensure that larger values correspond to higher collocability, the negative base-10 logarithm \eqn{-\log_{10} p}{-log10(p)} is returned, with negative sign iff \eqn{O_{11} < E_{11}}{O11 < E11}. Family-wise error rates are controlled with a Bonferroni correction to the p-values if \code{p.adjust=TRUE} (the default). Note that scores above 2 correspond to a significance level \eqn{p < .01}.
      }

      \item{\code{Fisher.pv}}{
        The \bold{p-values} of a one-sided \bold{Fisher's exact test} for independence of rows and columns in a contingency table, conditioned on the marginal frequencies. As usual for p-values, the negative base-10 logarithm \eqn{-\log_{10} p}{-log10(p)} is returned, which is always non-negative (keep in mind that the one-sided test only detects positive association). Family-wise error rates are controlled with a Bonferroni correction to the p-values if \code{p.adjust=TRUE} (the default).
      }

      \item{\code{simple.ll}}{
        Simple \bold{log-likelihood} (Evert 2008: 1225): \deqn{
          \pm 2 \left( O_{11}\cdot \log \frac{O_{11}}{E_{11}} - (O_{11} - E_{11}) \right)
        }{ ± 2 * ( O11 * log(O11 / E11) - (O11 - E11) ) }
        This measure provides a good approximation to the full log-likelihood measure (Evert 2008: 1235), but can be computed much more efficiently. Like \code{G2}, this measure computes \emph{signed} association scores and has a strong bias towards high co-occurrence frequency.
      }
      
      \item{\code{t}}{
        The \bold{t-score} association measure, which is popular for collocation identification in computational lexicography: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{O_{11}}}
        }{ (O11 - E11) / sqrt(O11) }
        T-score is known to filter out low-frequency data effectively.
      }
      
      \item{\code{X2}}{
        The \eqn{X^2} statistic of Pearson's \bold{chi-squared} test for independence of rows and columns in a contingency table, with Yates's correction applied: \deqn{
          \pm \frac{
            N \bigl( | O_{11}O_{22} - O_{12} O_{21} | - N/2 \bigr)^2
          }{
            R_1 R_2 C_1 C_2
          }
        }{ ± N * (|O12 * O22 - O12 * O21| - N/2)^2 / (R1 * R2 * C1 * C2) }
        This implementation computes \emph{signed} association scores, which are negative iff \eqn{O_{11} < E_{11}}{O11 < E11}. The formula above gives a more compact form of Yates's correction than the familiar sum over the four cells of the contingency table. See also \code{\link{chisq}}.
      }
      
      \item{\code{z}}{
        The \bold{z-score} association measure, based on a normal approximation to the binomial distribution of co-occurrence by chance: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{E_{11}}}
        }{ (O11 - E11) / sqrt(E11) }
        Z-score has a strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). Like pointwise MI, it should usually be combined with frequency thresholds on \eqn{f_1} and \eqn{f_2}.
      }
      
      \item{\code{Dice}}{
        The \bold{Dice coefficient} of association, which corresponds to the harmonic mean of the conditional probabilities \eqn{P(w_2|w_1)}{P(w2 | w1)} and \eqn{P(w_1|w_2)}{P(w1 | w2)}: \deqn{
          \frac{2 O_{11}}{R_1 + C_1}
        }{ 2 O11 / (R1 + C1) }
      }
      
      \item{\code{odds.ratio}}{
        Discounted log \bold{odds ratio}, an effect-size measure that is sensitive to small marginal frequencies: \deqn{
          \log \frac{(O_{11} + \frac12)(O_{22} + \frac12)}{(O_{12} + \frac12)(O_{21} + \frac12)}
        }{ log( (O11 + .5) * (O22 + .5) / (O12 + .5) * (O21 + .5) ) }
      }
      
      \item{\code{DP}}{
        The asymmetric \bold{Delta P} measure \eqn{\Delta P_{2|1} = P(w_2|w_1) - P(w_2|\neg w_1)}{Delta P_(2|1) = P(w2 | w1) - P(w2 | !w1)} proposed by Gries (2013: 143-144): \deqn{
          \frac{O_{11}}{R_1} - \frac{O_{21}}{R_2}
        }{ O11 / R1 - O21 / R2 }
      }
      
      \item{\code{LRC}}{
        The \bold{conservative LogRatio} (LRC) keyness measure (see \code{\link{keyness}}) can also be applied as an association measure. The implementation here always computes \code{PositiveLRC}, i.e. the lower boundary of a one-sided confidence interval for \eqn{\log_2 r}, in order to ensure a sensible and consistent scaling of the association score. Please keep in mind that negative scores do not necessarily indicate negative association (\eqn{O_{11} < E_{11}}{O11 < E11}), but rather lack of significant evidence for a positive association.
        
        Confidence intervals are determined at confidence level \code{conf.level} and adjusted with a Bonferroni correction if \code{p.adjust=TRUE} (the default). Note that LRC is applied to the rows of the contingency table rather than its columns, in order to obtain relative risk for \eqn{w_2} given \eqn{w_1}.
      }
    }
    
  }
  
  \subsection{User-defined association measures}{
  
     User-defined association measures can be applied by passing a suitable function in the \code{measure} argument. This function needs to be fully vectorised and will be applied to all co-occurrence pairs in the data set.
     
     It can use any of following arguments to access the contingency tables of observed and expected frequencies, following the notation of Evert (2008):
     \describe{
       \item{\code{O}, \code{E}}{observed and expected co-occurrence frequency}
       \item{\code{R1}, \code{R2}, \code{C1}, \code{C2}}{the row and column marginals of the contingency table}
       \item{\code{N}}{sample size}
       \item{\code{f}, \code{f1}, \code{f2}}{the frequency signature of each co-occurrence pair, a different notation for \eqn{f = O}{f = O}, \eqn{f_1 = R_1}{f1 = R1} and \eqn{f_2 = C_1}{f2 = C1}}
       \item{\code{O11}, \code{O12}, \code{O21}, \code{O22}}{the contingency table of observed frequencies}
       \item{\code{E11}, \code{E12}, \code{E21}, \code{E22}}{the contingency table of expected frequencies}
       \item{\code{param}}{a list with additional user-specified parameters, always including \code{conf.level} and \code{p.adjust}}
       \item{\code{\dots}}{must be specified to ignore unused arguments}
    }
    Except for \code{param}, all these arguments will be numeric vectors of the same length, and the function must return a numeric vector of the same length.
     
    For example, the built-in MI measure could also be implemented with the user function
    \preformatted{    my.MI <- function (O, E, ...) log2(O / E) }
    As a matter of fact, all built-in association measures are implemented in this way. The corresponding formulae can be obtained from the list \code{builtin.am}.

  }
  
  \subsection{Bonferroni correction}{
  
    If \code{p.adjust=TRUE}, statistical inference is corrected for multiple testing in order to control \bold{family-wise error rates}. This applies in particular to association measures based on confidence intervals (such as \code{LRC}) and versions that return (log) p-values (such as \code{G2-pv}). Note that the \code{G2} scores themselves are never adjusted.

    Family size \eqn{m} is automatically determined from the number of co-occurrence pairs processed in a single function call. Alternatively, the family size can be specified explicitly in the \code{p.adjust} argument, e.g. if a large data set is processed in multiple batches, or \code{p.adjust=FALSE} can be used to disable the correction.

    For the adjustment, a highly conservative Bonferroni correction \eqn{\alpha' = \alpha / m} is applied to significance levels. Since the large candidate sets and sample sizes often found in corpus linguistics tend to produce large numbers of false positives,  this conservative approach is considered to be useful.
  
  }

}

\seealso{

  \code{\link{prop.cint}}, which is used by some of the association measures based on confidence intervals; \code{\link{keyness}} for details on the LRC association measure (known as \code{PositiveLRC} there).

  For a gentle introduction to assocation measures, contingency tables, and different notions of co-occurrence
  see Evert (2008) or Unit 4 of the SIGIL course at \url{https://SIGIL.R-Forge.R-Project.org/}.
  The unit also includes a worked example carrying out several collocation analyses with real-life corpus data.

}
  
\examples{
## surface collocations with L2/R2 window
head(SurfaceColloc$cooc, 10) # table of word pairs with co-occurrence frequencies
head(SurfaceColloc$f1) # tables of marginal frequencies
head(SurfaceColloc$f2)

# add association scores (MI and p-values from log-likelihood test) to data frame
# keep in mind that we need to adjust expected frequencies for the total span size of 4 tokens
colloc <- transform(
  SurfaceColloc$cooc,
  MI = am.score(w1, w2, f, SurfaceColloc$f1, SurfaceColloc$f2, SurfaceColloc$N, 
                "MI", span.size=4),
  pval = am.score(w1, w2, f, SurfaceColloc$f1, SurfaceColloc$f2, SurfaceColloc$N, 
                  "G2.pv", span.size=4)
)
colloc <- colloc[order(colloc$MI, decreasing=TRUE), ] # ordered by MI scores
colloc

# apply significance filter
subset(colloc, pval > 3) # significant positive association at p < .001

# can also obtain labelled vector of scores
with(SurfaceColloc$cooc,
     am.score(w1, w2, f, SurfaceColloc$f1, SurfaceColloc$f2, SurfaceColloc$N, 
              "z", span.size=4, labels=TRUE))


## syntactic collocations from pre-computed frequency signatures (must specify named arguments!)
head(KrennPPV)
colloc <- transform(
  KrennPPV[, 1:9], # omit pre-computed association scores
  t = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="t"),
  OR = am.score(f=freq, f1=f.PP, f2=f.verb, N=N, measure="odds.ratio"))

all.equal(colloc$t, KrennPPV$t.score) # compare with pre-computed scores

colloc <- colloc[order(colloc$t, decreasing=TRUE), ]
head(colloc, 20) # top collocates (among all pairs)
subset(colloc, verb == "legen") # top collocates of verb "legen"


## collocations of "in charge of" with LRC as an association measure
colloc <- transform(
  BNCInChargeOf, # reconstruct frequency signature from O11, O12, C1, C2
  f = f.in,
  f1 = N.in,
  f2 = f.in + f.out,
  N = N.in + N.out
)
colloc <- transform(
  colloc,
  LRC = am.score("in charge of", collocate, f, f1, f2, N, "LRC"))
colloc <- colloc[order(colloc$LRC, decreasing=TRUE), ]
head(colloc, 20)
}

\references{

  \url{http://www.collocations.de/AM/}

  Evert, S. (2004). \emph{The Statistics of Word Cooccurrences: Word Pairs and Collocations}. Dissertation, Institut für maschinelle Sprachverarbeitung, University of Stuttgart. Published in 2005, URN urn:nbn:de:bsz:93-opus-23714. Available from \url{http://www.collocations.de/phd.html}.

  Evert, S. (2008). Corpora and collocations. In Lüdeling, A. and Kytö, M., editors, \emph{Corpus Linguistics. An International Handbook}, chapter 58, pages 1212–1248. Mouton de Gruyter, Berlin, New York. (\href{https://www.stephanie-evert.de/PUB/Evert2007HSK_extended_manuscript.pdf}{extended manuscript (PDF)})

  Evert, S. (2022). Measuring keyness. In \emph{Digital Humanities 2022: Conference Abstracts}, pages 202-205, Tokyo, Japan / online. 
  \url{https://osf.io/cy6mw/}
  
  Gries, S. T. (2013). 50-something years of work on collocations: What is or should be next …. \emph{International Journal of Corpus Linguistics}, 18(1):137–165.
 
}

\author{Stephanie Evert (\url{https://purl.org/stephanie.evert})}
\keyword{htest}