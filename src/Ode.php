<?php
declare(strict_types=1);


namespace iggyvolz\phpode;

use FFI;

class Ode
{
    private readonly FFI $ffi;

    public function __construct()
    {
        $this->ffi = FFI::cdef(file_get_contents(__DIR__ . "/cdef.h"), "libode.so");
        $this->ffi->dInitODE();
    }

    public function ffi(string $name, mixed ...$params): mixed
    {
        return $this->ffi->$name(...$params);
    }

    public function __destruct()
    {
        $this->ffi->dCloseODE();
    }
}