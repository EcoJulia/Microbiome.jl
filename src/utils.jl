macro forward_func(ex, fs)
    T, field = ex.args[1], ex.args[2].args[1]
    T = esc(T)
    fs = Base.Meta.isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end
