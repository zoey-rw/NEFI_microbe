#testing finding a scripts name. - this works
rstudioapi::getSourceEditorContext()$path

#just directory.
dirname(rstudioapi::getSourceEditorContext()$path)

#one level up.
