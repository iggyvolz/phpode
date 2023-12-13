<?php
declare(strict_types=1);


namespace iggyvolz\phpode;

use FFI\CData;

final class World
{
    private function __construct(public readonly Ode $ode, public readonly CData $cdata)
    {
    }

    public static function create(Ode $ode): self
    {
        return new self($ode, $ode->ffi("dWorldCreate"));
    }

    public function __destruct()
    {
        $this->ode->ffi("dWorldDestroy", $this->cdata);
    }

    public function setGravity(float $x, float $y, float $z): void
    {
        $this->ode->ffi("dWorldSetGravity", $this->cdata, $x, $y, $z);

    }

    public function step(float $len = 1): void
    {
        $this->ode->ffi("dWorldStep", $this->cdata, $len);
    }
}