

### copy mstoolkit from comet

```bash
cp -r ../Comet/MSToolkit/include .
find ../Comet/MSToolkit/src/expat-2.2.9/ -name "*.c" -exec cp {} src/expat-2.2.9/src \;
find ../Comet/MSToolkit/src/expat-2.2.9/ -name "*.h" -exec cp {} src/expat-2.2.9/include \;
find ../Comet/MSToolkit/src/zlib-1.2.11/ -name "*.c" -exec cp {} src/zlib-1.2.11/src \;
find ../Comet/MSToolkit/src/zlib-1.2.11/ -name "*.cc" -exec cp {} src/zlib-1.2.11/src \;
find ../Comet/MSToolkit/src/zlib-1.2.11/ -name "*.h" -exec cp {} src/zlib-1.2.11/include \;
cp -r ../Comet/MSToolkit/src/MSToolkit src
cp -r ../Comet/MSToolkit/src/mzIMLTools src
cp -r ../Comet/MSToolkit/src/mzParser src
```