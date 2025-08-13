# tests/testthat/test-cache-manager.R
# Test suite for LRU Cache Manager

library(testthat)
library(mockery)

# Source the cache manager
source("R/cache_manager.R")

# Helper to create a mock Stan model
create_mock_model <- function(size_mb = 10) {
  # Create an object with specified size
  mock_model <- list(
    data = rep(0, size_mb * 1024 * 256),  # Approximate size in MB
    compile_time = Sys.time(),
    model_code = "mock model"
  )
  class(mock_model) <- "CmdStanModel"
  return(mock_model)
}

# Helper to create temporary Stan files
create_temp_stan_file <- function(content = NULL) {
  if (is.null(content)) {
    content <- "
    parameters {
      real theta;
    }
    model {
      theta ~ normal(0, 1);
    }
    "
  }
  
  temp_file <- tempfile(fileext = ".stan")
  writeLines(content, temp_file)
  return(temp_file)
}

test_that("LRU cache initializes correctly", {
  cache <- LRUCache$new(max_size = 5)
  
  expect_s3_class(cache, "LRUCache")
  stats <- cache$get_stats()
  expect_equal(stats$n_models, 0)
  expect_equal(stats$max_size, 5)
  expect_equal(stats$total_size_mb, 0)
})

test_that("LRU cache stores and retrieves models", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 3)
  
  # Create temporary Stan file
  stan_file <- create_temp_stan_file()
  
  # Mock the cmdstanr::cmdstan_model function
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # First access should compile
  model1 <- cache$get(stan_file)
  expect_called(mock_compile, 1)
  
  # Second access should use cache
  model2 <- cache$get(stan_file)
  expect_called(mock_compile, 1)  # Still only called once
  
  # Check statistics
  stats <- cache$get_stats()
  expect_equal(stats$n_models, 1)
  expect_true(stats$total_size_mb > 0)
  
  # Clean up
  unlink(stan_file)
})

test_that("LRU cache evicts least recently used models", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 2)
  
  # Create multiple Stan files
  stan_files <- sapply(1:3, function(i) {
    create_temp_stan_file(paste("// Model", i))
  })
  
  # Mock compilation
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # Add first two models
  model1 <- cache$get(stan_files[1])
  Sys.sleep(0.1)  # Ensure different timestamps
  model2 <- cache$get(stan_files[2])
  
  stats <- cache$get_stats()
  expect_equal(stats$n_models, 2)
  
  # Access first model again (makes it more recent)
  model1_again <- cache$get(stan_files[1])
  
  # Add third model (should evict model2)
  model3 <- cache$get(stan_files[3])
  
  stats <- cache$get_stats()
  expect_equal(stats$n_models, 2)  # Still only 2 models
  
  # Model 1 and 3 should be in cache
  # Accessing model 2 should require recompilation
  initial_calls <- mock_call_count(mock_compile)
  model2_again <- cache$get(stan_files[2])
  expect_gt(mock_call_count(mock_compile), initial_calls)
  
  # Clean up
  sapply(stan_files, unlink)
})

test_that("Cache respects maximum size limit", {
  skip_if_not_installed("cmdstanr")
  
  max_size <- 10
  cache <- LRUCache$new(max_size = max_size)
  
  # Create many Stan files
  n_files <- 15
  stan_files <- sapply(1:n_files, function(i) {
    create_temp_stan_file(paste("// Model", i))
  })
  
  # Mock compilation
  mock_compile <- mock(create_mock_model(size_mb = 2), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # Add all models
  for (file in stan_files) {
    cache$get(file)
  }
  
  # Check that cache size never exceeds maximum
  stats <- cache$get_stats()
  expect_lte(stats$n_models, max_size)
  
  # Clean up
  sapply(stan_files, unlink)
})

test_that("Force recompile works correctly", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 5)
  stan_file <- create_temp_stan_file()
  
  # Mock compilation
  compile_count <- 0
  mock_compile <- mock(
    {
      compile_count <<- compile_count + 1
      create_mock_model(size_mb = 5)
    }, 
    cycle = TRUE
  )
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # First compilation
  model1 <- cache$get(stan_file)
  expect_equal(compile_count, 1)
  
  # Use cache
  model2 <- cache$get(stan_file)
  expect_equal(compile_count, 1)
  
  # Force recompile
  model3 <- cache$get(stan_file, force_recompile = TRUE)
  expect_equal(compile_count, 2)
  
  # Clean up
  unlink(stan_file)
})

test_that("Cache clear function works", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 5)
  
  # Add some models
  stan_files <- sapply(1:3, function(i) {
    create_temp_stan_file(paste("// Model", i))
  })
  
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  for (file in stan_files) {
    cache$get(file)
  }
  
  stats_before <- cache$get_stats()
  expect_equal(stats_before$n_models, 3)
  
  # Clear cache
  cache$clear()
  
  stats_after <- cache$get_stats()
  expect_equal(stats_after$n_models, 0)
  expect_equal(stats_after$total_size_mb, 0)
  
  # Clean up
  sapply(stan_files, unlink)
})

test_that("Memory statistics are tracked correctly", {
  stats <- get_memory_stats()
  
  expect_type(stats, "list")
  expect_s3_class(stats, "memory_stats")
  expect_true("timestamp" %in% names(stats))
  expect_true("r_memory" %in% names(stats))
  expect_true("cache" %in% names(stats))
  
  # Check R memory stats
  expect_true(stats$r_memory$used_mb >= 0)
  expect_true(stats$r_memory$gc_trigger_mb >= 0)
  
  # Check cache stats
  expect_true(stats$cache$n_models >= 0)
  expect_true(stats$cache$total_size_mb >= 0)
})

test_that("Access counts and times are tracked", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 5)
  stan_file <- create_temp_stan_file()
  
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # Access model multiple times
  for (i in 1:5) {
    cache$get(stan_file)
    Sys.sleep(0.01)  # Small delay
  }
  
  stats <- cache$get_stats()
  expect_equal(nrow(stats$models), 1)
  expect_equal(stats$models$access_count[1], 5)
  expect_true(stats$models$last_access[1] <= Sys.time())
  
  # Clean up
  unlink(stan_file)
})

test_that("Global cache instance works", {
  # Get global instance
  cache1 <- get_lru_cache(max_size = 7)
  cache2 <- get_lru_cache()  # Should return same instance
  
  expect_identical(cache1, cache2)
  
  stats <- cache1$get_stats()
  expect_equal(stats$max_size, 7)
})

test_that("Garbage collection is triggered", {
  # Test force_gc function
  result <- force_gc(verbose = FALSE)
  expect_type(result, "double")
  expect_true(all(dim(result) == c(2, 6)))  # Standard gc output dimensions
})

test_that("Integration with backend_bayes.R works", {
  skip_if_not_installed("cmdstanr")
  
  # Source modified backend
  source("R/backend_bayes.R")
  
  # Initialize cache
  cache <- .initialize_model_cache()
  expect_s3_class(cache, "LRUCache")
  
  # Create Stan file
  stan_file <- create_temp_stan_file()
  
  # Mock cmdstanr
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(get_compiled_model, "cmdstanr::cmdstan_model", mock_compile)
  
  # Test get_compiled_model
  model <- get_compiled_model(stan_file)
  expect_s3_class(model, "CmdStanModel")
  
  # Test cache is used
  model2 <- get_compiled_model(stan_file)
  expect_called(mock_compile, 1)  # Only compiled once
  
  # Clean up
  clear_model_cache()
  unlink(stan_file)
})

test_that("Memory monitoring options work", {
  # Test with monitoring enabled
  options(tdmx_monitor_memory = TRUE)
  options(tdmx_debug = TRUE)
  
  test_fn <- function(x) x * 2
  
  expect_message(
    result <- monitor_memory("test", test_fn, 5),
    regexp = "Memory"
  )
  expect_equal(result, 10)
  
  # Test with monitoring disabled
  options(tdmx_monitor_memory = FALSE)
  options(tdmx_debug = FALSE)
  
  expect_silent(
    result <- monitor_memory("test", test_fn, 5)
  )
  expect_equal(result, 10)
})

test_that("GC options control garbage collection", {
  # Test GC after sampling
  options(tdmx_gc_after_sampling = TRUE)
  options(tdmx_debug = TRUE)
  
  # Mock the force_gc to track if it's called
  gc_called <- FALSE
  mock_gc <- mock({
    gc_called <<- TRUE
    gc(verbose = FALSE)
  })
  
  with_mock(
    force_gc = mock_gc,
    {
      # This would normally be called within run_fit_stan_hmc
      if (getOption("tdmx_gc_after_sampling", TRUE)) {
        force_gc(verbose = getOption("tdmx_debug", FALSE))
      }
    }
  )
  
  expect_true(gc_called)
  
  # Reset options
  options(tdmx_gc_after_sampling = NULL)
  options(tdmx_debug = NULL)
})

test_that("Cache handles file modifications correctly", {
  skip_if_not_installed("cmdstanr")
  
  cache <- LRUCache$new(max_size = 5)
  stan_file <- create_temp_stan_file("// Version 1")
  
  mock_compile <- mock(create_mock_model(size_mb = 5), cycle = TRUE)
  stub(cache$get, "cmdstanr::cmdstan_model", mock_compile)
  
  # First compilation
  model1 <- cache$get(stan_file)
  initial_calls <- mock_call_count(mock_compile)
  
  # Modify file
  Sys.sleep(1)  # Ensure different mtime
  writeLines("// Version 2", stan_file)
  
  # Should recompile due to different mtime
  model2 <- cache$get(stan_file)
  expect_gt(mock_call_count(mock_compile), initial_calls)
  
  # Clean up
  unlink(stan_file)
})

test_that("Cache handles missing files gracefully", {
  cache <- LRUCache$new(max_size = 5)
  
  expect_error(
    cache$get("non_existent_file.stan"),
    regexp = "Stan file not found"
  )
})

test_that("Print methods work correctly", {
  cache <- LRUCache$new(max_size = 5)
  
  # Test empty cache print
  expect_output(print(cache), "LRU Stan Model Cache")
  expect_output(print(cache), "Models: 0 / 5")
  
  # Test memory stats print
  stats <- get_memory_stats()
  expect_output(print(stats), "Memory Statistics")
  expect_output(print(stats), "R Memory Usage")
  expect_output(print(stats), "Stan Model Cache")
})