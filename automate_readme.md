## automate_readme.md

Set up env vars automatically:

1. install direnv, 
2. add hook to .bashrc file:
>`eval "$(direnv hook bash)"`

[Helpful blogpost](https://kellner.io/direnv.html)

## Quick direnv Demo

>Show that the FOO environment variable is not loaded.
```
$ echo ${FOO-nope}
nope
```

>Create a new .envrc. This file is bash code that is going to be loaded by direnv.
```
$ echo export FOO=foo > .envrc
.envrc is not allowed
```

>The security mechanism didn't allow to load the .envrc. Since we trust it, let's allow its execution.
```
$ direnv allow .
direnv: reloading
direnv: loading .envrc
direnv export: +FOO
```

>Show that the FOO environment variable is loaded.
```
$ echo ${FOO-nope}
foo
```

>Exit the project
```
$ cd ..
direnv: unloading
```

>And now FOO is unset again
```
$ echo ${FOO-nope}
nope
```