This is a TypeScript version of some interpolation Java classes of the
Apache Commons Math library. It currently contains the following main
classes and related classes and interfaces:

 - AkimaSplineInterpolator
 - SplineInterpolator
 - LinearInterpolator

The source code has been manually translated from Java to Typescript.
Unneeded class methods have been omitted.

Copyright 2016 Christian d'Heureuse, Inventec Informatik AG, Zurich, Switzerland
www.source-code.biz, www.inventec.ch/chdh

-----------------------------------------------------------------------------

Copyright 2001-2016 The Apache Software Foundation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

-----------------------------------------------------------------------------
* to build I used `npm install && node_modules/.bin/tsc --declaration ApacheCommonsMathInterpolation.ts && node_modules/.bin/uglifyjs ApacheCommonsMathInterpolation.js >ApacheCommonsMathInterpolation.min.js`
