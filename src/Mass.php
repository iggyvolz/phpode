<?php
declare(strict_types=1);


namespace iggyvolz\phpode;

use FFI;
use FFI\CData;

class Mass
{
    private function __construct(private Ode $ode, public readonly CData $cdata)
    {
    }

    public static function create(Ode $ode): self
    {
        return new self($ode, $ode->ffi("new", "dMass"));
    }

    public static function box(Ode $ode, float $density, float $lx, float $ly, float $lz): self
    {
        $self = self::create($ode);
        $self->setBox($density, $lx, $ly, $lz);
        return $self;
    }

    public function setBox(float $density, float $lx, float $ly, float $lz): void
    {
        $this->ode->ffi("dMassSetBox", FFI::addr($this->cdata), $density, $lx, $ly, $lz);
    }
}