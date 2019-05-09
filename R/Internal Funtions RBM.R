# Internal Functions from https://github.com/TimoMatzen/RBM

# Restricted Boltzmann Machine#
RBM <- function (x, y, n.iter = 100, n.hidden = 30, learning.rate = 0.1,
                 plot = FALSE, size.minibatch = 10, momentum = 0.5, lambda = 0.001) {
  # Trains a Restricted Boltzmann Machine on binary data, either supervised
  # or unsupervised. This function uses contrastive diversion with k = 1 for training the system.
  #
  # Args:
  #   x: A matrix with binary features of shape samples * features.
  #   y: A matrix with labels for the data, only when training a classification RBM. (Optional)
  #   n.iter: Defines the number of epochs to run contrastive diversion.
  #   n.hidden: The number of nodes in the hidden layer.
  #   learning.rate: The learning rate, alpha, for training the system.
  #   plot: Whether to plot the learning progress of the weights
  #   size.minibatch: The size of the minibatches used for training.
  #   momentum: Speeds up the gradient descent learning.
  #   lambda: The sparsity penalty lambda to prevent the system from overfitting.
  #
  # Returns:
  #   A list with the trained weights of the RBM that can be used for the predict RBM function when supervised learning was applied
  #   or the ReconstructRBM function to reconstruct data with the model.

  # Check whether n.iter is devicable by ten and if so initialize plot.epoch:
  if (plot == TRUE) {
    if ((n.iter %% 10) == 0) {
      # Plot at each n.iter/10
      plot.epoch <- n.iter/10
    } else {
      # Warn user and turn plotting off
      print ('Number of iterations was not devicable by ten: plots are turned off')
      plot <- FALSE
      plot.epoch <- 0
    }
  } else {
    # Set plot.epoch to FALSE
    plot.epoch <- FALSE
  }

  # Some checks
  if (!is.matrix(x)) {
    print('Data was not in a matrix, converted data to a matrix')
    x <- as.matrix(x)
  }
  if (any(!is.numeric(x))) {
    stop('Sorry the data has non-numeric values, the function is terminated')
  }
  if (n.iter > 10000) {
    print("Number of epochs for > 10000, could take a while to fit")
  }
  if (!missing(y)) {
    if (any(!is.numeric(y))) {
      stop('Sorry the labels have non-numeric values, the function is terminated')
    }
    if (any(!is.finite(y))) {
      stop('Sorry this function cannot handle NAs or non-finite label values')
    }
    if (length(y) != nrow(x)) {
      stop('Labels and data should be equal for supervised RBM: try training an unsupervised RBM')
    }
  }
  if (any(!is.finite(x))) {
    stop('Sorry this function cannot handle NAs or non-finite data')
  }
  if (size.minibatch > 100) {
    print('Sorry the size of the minibatch is too long: resetting to 10')
    size.minibatch <- 10
  }
  if (size.minibatch > 20) {
    print('Large minibatch size, could take a long time to fit model')
  }
  if (min(x) < 0 | max(x) > 1) {
    stop('Sorry the data is out of bounds, should be between 0 and 1')
  }
  if( length(dim(x)) < 2 ) {
    stop("Dimensions of the data were not right, should be of shape n.features * n.samples")
  }
  if(ncol(x) > nrow(x)) {
    print('Less data than features, this will probably result in a bad model fit')
  }

  # Initialize the weights, n.features * n.hidden with values from gaussian distribution
  weights <- matrix(rnorm(ncol(x) * n.hidden, 0, .01), nrow = ncol(x), ncol = n.hidden)
  # Initialize the momentum_speed matrix
  momentum_speed_x <- matrix(0, nrow = ncol(x) + 1, ncol = n.hidden + 1)

  # Add bias to weights
  weights <- cbind(0, weights)
  weights <- rbind(0, weights)

  # Add 1 for the bias to x
  x <- cbind(1, x)

  # Initialize the labels, weights and bias for the labels if supervised = TRUE
  if (!missing(y)) {
    # Get all the unique labels in y
    labels <- unique(y)
    # Get the indexes of each unique label in y
    idx <- vector('list', length = length(labels))
    # Save indexes
    for (i in 1:length(labels)) {
      idx[[i]]<- which(y == labels[i])
    }
    # Make binarized vectors of the labels
    y <- LabelBinarizer(y)
    # Add one term for the bias
    y <- cbind(1, y)

    # Create the y weights matrix
    y.weights <- matrix(rnorm(length(labels) * n.hidden, 0, 01), nrow = length(labels), ncol = n.hidden)
    # Add momentum speed matrix
    momentum_speed_y <- matrix(0, nrow = length(labels) + 1, ncol = n.hidden + 1)

    # add bias to weights
    y.weights <- cbind(0, y.weights)
    y.weights <- rbind(0, y.weights)

  }
  # PLot the untrained weights
  if(plot == TRUE){
    # Set plotting margins
    par(mfrow = c(3,10), mar = c(3,1,1,1))
    plot.weights <- weights[-1, -1]
    if (n.hidden > 30) {
      # Warn user that only a sample of the hidden nodes will plotted
      print('n.hidden > 30, only plotting a sample of the invisible nodes')
      # Take sample
      samp.plot <- sample(1:n.hidden, 30)
      # Remove weights for plotting
      for(i in samp.plot) {
        # Plot weights
        image(matrix(plot.weights[, i], nrow = sqrt(ncol(x)-1)), col=grey.colors(255))
        title(main = paste0('Hidden node ', i), font.main = 4)
        # Initialize counter for the plotting:
        plot.counter <- 0
      }
    } else {
      for(i in 1:n.hidden) {
        # Plot weights
        image(matrix(plot.weights[, i], nrow = sqrt(ncol(x)-1)), col=grey.colors(255))
        title(main = paste0('Hidden node ', i), font.main = 4)
        # Initialize counter for the plotting:
        plot.counter <- 0
      }
    }
  }
  plot.counter <- 0
  # Start contrastive divergence, k = 1
  for (i in 1:n.iter){
    if (missing(y)) {
      # Sample minibatch from x, unsupervised
      samp <- sample(1:nrow(x), size.minibatch, replace = TRUE)
    } else {
      # Pick balanced labels
      samp <- rep(0,size.minibatch)
      p <- 1
      for (i in 1 : size.minibatch){
        samp[p]<- sample(idx[[p]], 1)
        p <- p + 1
        if (p == length(labels) +1) {
          # Reset counter
          p <- 1
        }
      }
    }
    plot.counter <- plot.counter + 1
    # At iteration set visible layer to random sample of train:
    V0 <- x[samp, ,drop = FALSE]
    if (missing(y)) {
      # Calculate gradients
      grads <- CD(V0, weights)
    } else {
      # Calculate gradients
      grads <- CD(V0, weights, y[samp,,drop = FALSE], y.weights)
    }
    # Update the momentum speed
    momentum_speed_x <- momentum * momentum_speed_x + ((grads$grad.weights - (lambda * weights))/ size.minibatch)

    # Update weights and bias
    weights <- weights + (learning.rate * momentum_speed_x)

    if (!missing(y)) {
      # Update momentum speed
      momentum_speed_y <- momentum * momentum_speed_y + ((grads$grad.y.weights - (lambda * y.weights))/ size.minibatch)


      # Update weights and bias
      y.weights <- y.weights + (learning.rate * momentum_speed_y)
    }
    # Plot learning of hidden nodes at every plot.epoch:
    if(plot.counter == plot.epoch & plot == TRUE) {
      # Create margins
      par(mfrow = c(3,10), mar = c(3,1,1,1))
      # Remove bias for plottingun
      plot.weights <- weights[-1, -1]
      if (n.hidden > 30) {
        for(i in samp.plot) {
          image(matrix(plot.weights[, i], nrow = sqrt(ncol(x)-1)), col=grey.colors(255))
          title(main = paste0('Hidden node ', i), font.main = 4)
        }
      } else {
        for(i in 1:n.hidden) {
          image(matrix(plot.weights[, i], nrow = sqrt(ncol(x)-1)), col=grey.colors(255))
          title(main = paste0('Hidden node ', i), font.main = 4)
        }
      }
      # Reset the plot counter:
      plot.counter <- 0
    }
  }
  # Return list with the matrices of trained weights and bias terms
  if (!missing(y)) {
    return(list('trained.weights' = weights,'trained.y.weights' = y.weights))
  } else {
    return(list('trained.weights' = weights))
  }
}


# Function to calculate hidden layer from data
#
# @keyword internal
#
# Function for calculating hidden layer:
VisToHid <- function(vis, weights, y, y.weights) {
  # Function for calculating a hidden layer.
  #
  # Args:
  #   vis: Visual layer, or hidden layer from previous layer in DBN
  #   weights: Trained weights including the bias terms (use RBM)
  #   y: Label vector if only when training an RBM for classification
  #   y.weights: Label weights and bias matrix, only neccessary when training a RBM for classification
  #
  # Returns:
  #   Returns a hidden layer calculated with the trained RBM weights and bias terms.
  #
  # Initialize the visual, or i-1 layer
  V0 <- vis
  if ( is.null(dim(V0))) {
    # If visual is a vector create matrix
    V0 <- matrix(V0, nrow= length(V0))
  }
  if(missing(y) & missing(y.weights)) {
    # Calculate the hidden layer with the trained weights and bias
    H <- 1/(1 + exp(-( V0 %*% weights)))
  } else {
    Y0 <- y
    H <- 1/(1 + exp(- ( V0 %*% weights + Y0 %*% y.weights)))
  }
  return(H)
}

# Function for reconstructing data from a hidden layer
#
# @keyword internal
# Function for reconstructing visible layer:
HidToVis <- function(inv, weights, y.weights) {
  # Function for reconstructing a visible layer.
  #
  # Args:
  #   inv: Invisible layer
  #   vis.bias: Trained visible layer bias (use RBM)
  #   weights: Trained weights (use RBM)
  #   y.weights: Label weights, only nessecessary when training a classification RBM.
  #
  # Returns:
  #   Returns a vector with reconstructed visible layer or reconstructed labels.
  #
  if(missing(y.weights)) {
    # Reconstruct only the visible layer when y.weights is missing
    V <- 1/(1 + exp(-(   inv %*% t(weights)) ))
    return(V)
  } else {
    # Reconstruct visible and labels if y.weights
    Y <- 1/(1 + exp(-( inv %*% t(y.weights))))
    return(Y)
  }
}

# Function for doing contrastive divergence CD
#
# @keyword internal
CD <- function(vis, weights, y, y.weights) {
  # Function for doing k=1 contrastive divergence
  #
  # Args:
  #   vis: visible layer values vector of shape n_features * 1
  #   weights: weights vector of shape n_features * n_hidden
  #   vis.bias: bias of the visible layer
  #   inv.bias: bias of the invisible layer
  #   y: labels, only used when provided
  #   y.weigths: label weights of shape n_labels * n_hidden, only used when provided
  #   y.bias: bias term for the labels of shape n_features * 1, only used when provided
  #
  # Returns:
  #   A list with all gradients for the bias and weights; adds label bias and weights if y is provided
  #
  # Start positive phase
  if (missing(y) & missing(y.weights)) {
    # Calculate hidden layer
    H0 <- VisToHid(vis, weights)
    H0[,1] <- 1
  } else {
    # Add a layer with labels if y is provided
    H0 <- VisToHid(vis, weights, y, y.weights)
    H0[,1] <- 1
  }
  # Binarize the hidden layer:
  unif  <- runif(nrow(H0) * (ncol(H0)))
  H0.states <- H0 > matrix(unif, nrow=nrow(H0), ncol= ncol(H0))

  # Calculate positive phase, we always use the probabilities for this
  pos.phase <- t(vis) %*% H0
  if (!missing(y)) {
    pos.phase.y <- t(y) %*% H0
  }
  # Start negative  phase
  # Reconstruct visible layer
  V1 <- HidToVis(H0.states, weights)
  # Set the bias unit to 1
  V1[,1] <- 1

  if (missing(y) & missing(y.weights) ) {
    # Reconstruct hidden layer unsupervised, no need to fix the bias anymore
    H1 <- VisToHid(V1, weights)
  } else {
    # Reconstruct labels if y is provided
    Y1 <- HidToVis(H0, weights,  y.weights )
    # Set the bias unit to 1
    Y1[,1] <- 1
    # Reconstruct hidden layer supervised, no need to fix the bias anymore
    H1 <- VisToHid(V1, weights, Y1, y.weights)
  }
  # Calculate negative associations, we alway use the probabilities for this:
  neg.phase <- t(V1) %*% H1
  if (!missing(y) & !missing(y.weights)) {
    # Calculate negative phase y
    neg.phase.y <- t(Y1) %*% H1
  }
  ## Calculate the gradients
  # Calculate gradients for the weights:
  grad.weights <- pos.phase - neg.phase

  if (!missing(y) & !missing(y.weights)) {
    # Calculate gradients for y.weigths
    grad.y.weights <- pos.phase.y - neg.phase.y

    # Return list with  gradients supervised
    return(list('grad.weights' = grad.weights,'grad.y.weights' = grad.y.weights))
  } else {
    # Return list with gradients unsupervised
    return(list('grad.weights' = grad.weights  ))
  }
}

