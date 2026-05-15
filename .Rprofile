local({
  # Prefer reticulate's managed uv environment unless the user explicitly
  # pinned a Python interpreter for this session.
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
    Sys.setenv(RETICULATE_PYTHON = "managed")
  }
})
