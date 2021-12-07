# Based on Lab 2:
# Exponential Random Graph Models (ERGMs)
#
# Install packages below if you do not have them:
# -------------------------------------------------
# install.packages("statnet")
# install.packages("igraph")
# install.packages("texreg")
detach(package:igraph)
library(statnet)
?statnet
# -------------------------------------------------------------------------------------------------
# Set the working directory
# Session > Set Working Directory > To Source File Location
# -------------------------------------------------------------------------------------------------
list.files() # List the files in the current working directory to see if you're in the right directory

# ----------------------------------------------------------------------------------------------------
######################## PART I: Building and Visualizing the Networks ########################
# ----------------------------------------------------------------------------------------------------
# Dependent variable:
# Responses to the question:
# " List up to 5 employees at <your company> who are the most capable of getting buy-in from you and changing your behaviors and opinions."
# Note that participants were limited to selecting at most five respondents.
link <- read.csv("link.csv")
# View the first rows of the edgelist to make sure it imported correctly:
head(link)
# Convert the edgelist to a network object in statnet format:
replies <- as.network.matrix(link, matrix.type = "edgelist") 
replies # View a summary of the network object

# Independent variables:
# Load node attributes, and store them in the buyIn network object we have created
set.vertex.attribute(replies, "is_pos", read.csv("is_pos.csv")$is_pos)
set.vertex.attribute(replies, "is_neg", read.csv("is_neg.csv")$is_neg)
set.vertex.attribute(replies, "is_extreme", read.csv("is_extreme.csv")$is_extreme)
set.vertex.attribute(replies, "pos", read.csv("pos.csv")$pos)
set.vertex.attribute(replies, "neg", read.csv("neg.csv")$neg)
set.vertex.attribute(replies, "neu", read.csv("neu.csv")$neu)
set.vertex.attribute(replies, "compound", read.csv("compound.csv")$compound)
set.vertex.attribute(replies, "length", read.csv("length.csv")$length/100)
set.vertex.attribute(replies, "upvote", read.csv("upvote.csv")$upvote/1000)
set.vertex.attribute(replies, "upvote_thres", read.csv("upvote_thres.csv")$upvote_thres)

replies # These five variables should now be listed as vertex attributes when viewing the summary of the network
# Check the values for all of the node variables
get.vertex.attribute(replies,"pos")
get.vertex.attribute(replies,"neg")
get.vertex.attribute(replies,"neu")
get.vertex.attribute(replies,"compound")
get.vertex.attribute(replies,"length")

get.vertex.attribute(replies,"upvote")

# ---------------------------------------------------------------------------------------
# Basic information
# ---------------------------------------------------------------------------------------
summary(replies)                              
network.size(replies)                        
betweenness(replies)                         
isolates(replies)                            

# ---------------------------------------------------------------------------------------
# Visualize networks
# ---------------------------------------------------------------------------------------
library('igraph')

# Set default plot options
igraph_options(vertex.size = 5, vertex.color = 'grey', # vertex.size changes the size of nodes; vertex.color changes the color of nodes
               edge.color='gray80', edge.arrow.size=.4, # edge.color changes the color of ties; edge.arrow.size changes the size of tie arrow heads
               vertex.label = NA)                       # vertex.label = NA specifies not to display vertex labels in the plot


replies_igraph <- graph.adjacency(as.matrix.network(replies))
replies_igraph <- set_vertex_attr(replies_igraph,"is_pos",value = read.csv("is_pos.csv")$is_pos)
replies_igraph <- set_vertex_attr(replies_igraph,"upvote_thres",value = read.csv("upvote_thres.csv")$upvote_thres)
net_layout <- layout_with_fr(replies_igraph)
plot(replies_igraph, layout=net_layout, edge.color='black', vertex.label = NA)


V(replies_igraph)$color = ifelse (V(replies_igraph)$is_pos ==1, " orange ", " grey ")
plot(replies_igraph, layout=net_layout, edge.color='black', vertex.color = ifelse (V(replies_igraph)$upvote_thres ==1, " red ", " grey ") )


# -------------------------------------------------------------------------------------------------
######################## PART II: Build the ERGM models ########################
#
# R vignette for more details: https://cran.r-project.org/web/packages/ergm/ergm.pdf
# -------------------------------------------------------------------------------------------------
# Remove the 'igraph' package from your environment. 
detach(package:igraph)
library(statnet)
options(ergm.loglik.warn_dyads=FALSE) #Whether or not a warning should be issued when sample space constraints render the observed number of dyads ill-defined

# Ergm Terms are statistics: They are some deterministic function of the ties, node attributes, and edge covariates of a network.
help('ergm-terms') # Documentation that contains definitions for all of the terms we are using
                   # ex. what does "mutual" test and how is it calculated
                   # You will want to click the link for "Terms used in Exponetnail Family Random Graph Models" if offered it
# We will use the ergm-terms to perform hypothesis testing using ERGMs
# But we can note that any of the ERGM terms can also be examined directly for your obsved network, by creating a formula in R
summary (replies ~ edges + mutual) # Number of edges, number of reciprocated dyads
summary (replies ~ odegree(0:5))   # Outdegree distribution. Remember, respondents could nominate at most five employees
summary (replies ~ idegree(0:65))  # Indegree distribution.
# This type of analysis can be helpful for understanding your network, as well as troubleshooting issues with ERGM regression


model2 <- ergm(replies ~  # This model will be slower to estimate than model 1
                        # Expect roughly 2-7 minutes. If it gets stuck for longer than that, try hitting "stop" and re-running it
               # Structural patterns
               # edges
               #mutual
               #gwidegree(log(2), fixed = T)                 # Inverted preferential attachment (indegree)
               #gwodegree(2, fixed = T, cutoff = 100)              # Inverted preferential attachment (outdegree)
               #dgwesp(log(2), type = "OTP", fixed = T, cutoff =5)    # A modified version of Outgoing Two Path(i->j + i->k->j) structures. Geometrically weighted version of transitivity
               # Node attribute effects
               + nodematch("is_pos")                                   # Homophily on a categorical variable 
               + nodematch("upvote_thres")
               #+ nodematch("is_extreme") 
               #+ nodemix("is_pos", base = c(0, 1))
               #+ nodematch("is_neg")
               #+ nodematch("pos")                           # Mixing matrix of all different combinations of node attributes (ex. A -> A ties, A-> B ties, B -> A ties, B -> B ties). 
               #+ nodematch("neg") 
               #+ nodeicov("office")                                    # Covariance between in-degree of nodes and attributes of nodes
               + nodeocov("pos")                                    # Covariance between out-degree of nodes and attributes of nodes
               + nodeocov("neg") 
               + nodeocov("neu")
               + nodeocov("length")
               #+ nodeocov("upvote")
               + nodeocov("compound")
               + nodeocov("is_extreme")
               + nodeocov("upvote_thres")
               , constraints =~ bd(maxout=100)                           # This constraint enforces the maximum outdegree is 5
               # Control settings for MCMC-MLE algorithm
               , control=snctrl(MCMLE.maxit=50, MCMC.interval = 100)
               , verbose=T
) 
summary(model2) 

mcmc.diagnostics(model2)

# Easy side-by-side model comparison:
library(texreg)
# -------------------------------------------------------------------------------------------------
######################## PART III: Model diagnostics ########################
# 
# Change 'file_name' and 'model' variables to refer to the desired model (model1, or model2)
# -------------------------------------------------------------------------------------------------
file_name  <- 'model_diagnostics_1.pdf'
model      <- model2           # Insert one of the models here! (model1 or model2)
                               # If you change / re-estimate your models, make sure you re-run the above command
                               # before you run any of the below commands

pdf(file_name)  # Saves the model diagnostic as a PDF - look for this in your current working directory
mcmc.diagnostics(model)         # Performs the markov chain monte carlo diagnostics
dev.off()                       # Closes the pdf saving process.

# -------------------------------------------------------------------------------------------------
# Goodness of fit test
# Check how well the estimated model captures certain features of the observed network, for example triangles in the network.
# -------------------------------------------------------------------------------------------------
# This first command simulates 100 networks.
# These networks, if we use sufficient burnin steps in the markov chain used to generate them,
# may be thought of as random samples from the joint probability distribution that is our fitted ERGM.
sim <- simulate(model2, burnin=100000, interval=100000, nsim=100, verbose=T)  # Uses the ergm model to simulate a null model

# Plot the first of the simulated networks
sim1_net <- igraph::graph.adjacency(as.matrix.network(sim[[1]]))
igraph::plot.igraph(sim1_net,layout=net_layout,edge.color="brown",  
                    vertex.color = 'grey',edge.arrow.size=.4)                                                               

# Plot the 10th simulated network
sim10_net <- igraph::graph.adjacency(as.matrix.network(sim[[10]]))
igraph::plot.igraph(sim10_net,layout=net_layout,edge.color="purple",  
                    vertex.color = 'grey',edge.arrow.size=.4)                                                                 

# -------------------------------------------------------------------------------------------------
# Extract the number of triangles from each of the 100 samples and
# compare the distribution of triangles in the sampled networks with the observed network
# -------------------------------------------------------------------------------------------------
model2.tridist <- sapply(1:100, function(x) summary(sim[[x]] ~triangle)) # Extracts the tiangle data from the simulated networks
hist(model2.tridist,xlim=c(0,1000),breaks=10)                             # Plots that triangle distribution as a histogram, change xlim to change the x-axis range if necessary
replies.tri <- summary(replies ~ triangle)                                  # Saves the CRIeq triangle data from the summary to the CRI.eq variable
replies.tri
arrows(replies.tri,20, replies.tri, 5, col="red", lwd=3)                    # Adds an arrow to the plotted histogram
c(obs=replies.tri,mean=mean(model2.tridist),sd=sd(model2.tridist),
  tstat=abs(mean(model2.tridist)-replies.tri)/sd(model2.tridist))

# -------------------------------------------------------------------------------------------------
# Test the goodness of fit of the model
# Compiles statistics for these simulations as well as the observed network, and calculates p-values 
# -------------------------------------------------------------------------------------------------
# This first command runs goodness of fit testing
# It may take a second for this command to run.
gof <- gof(model2, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
# If you run below and then wouldn't see the plot, trypar(mar=c(2,2,2,2))
plot(gof)           # Plot the goodness of fit
                    # Note: this should produce five separate plots that you should look through.
gof                 # Display the goodness of fit info in the console

