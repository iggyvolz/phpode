<?php
declare(strict_types=1);


namespace iggyvolz\phpode;

final readonly class Vec3
{
    public function __construct(public float $x, public float $y, public float $z)
    {
    }

    public function __toString(): string
    {
        return "($this->x, $this->y, $this->z)";
    }
}