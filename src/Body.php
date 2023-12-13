<?php
declare(strict_types=1);


namespace iggyvolz\phpode;

use FFI;
use FFI\CData;

final class Body
{
    private Ode $ode;

    private function __construct(private World $world, private CData $cdata)
    {
        $this->ode = $this->world->ode;
    }

    public static function create(World $world): self
    {
        return new self($world, $world->ode->ffi("dBodyCreate", $world->cdata));
    }

    public function setMass(Mass $mass): void
    {
        $this->ode->ffi("dBodySetMass", $this->cdata, FFI::addr($mass->cdata));
    }

    public function __destruct()
    {
        $this->ode->ffi("dBodyDestroy", $this->cdata);
    }

    public function getPosition(): Vec3
    {
        $ptr = $this->ode->ffi("dBodyGetPosition", $this->cdata);
        return new Vec3($ptr[0], $ptr[1], $ptr[2]);
    }
}