<?php

use iggyvolz\phpode\Body;
use iggyvolz\phpode\Mass;
use iggyvolz\phpode\Ode;
use iggyvolz\phpode\World;

require_once __DIR__ . "/vendor/autoload.php";
$ode = new Ode();
$world = World::create($ode);
$world->setGravity(0, 0, -9.81);
$body = Body::create($world);
$body->setMass(Mass::box($ode, 1, 1, 1, 1));
echo "Initial box position: " . $body->getPosition() . PHP_EOL;
for ($i = 0; $i < 5; $i++) {
    $world->step();
    echo "Box position after step $i: " . $body->getPosition() . PHP_EOL;
}
