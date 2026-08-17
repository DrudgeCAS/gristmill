# Changelog

## [0.9.0](https://github.com/DrudgeCAS/gristmill/compare/gristmill-v0.8.0...gristmill-v0.9.0) (2026-08-17)


### ⚠ BREAKING CHANGES

* Migrate to scikit build system ([#22](https://github.com/DrudgeCAS/gristmill/issues/22))

### Bug Fixes

* Apply the three review fixes to the search ([5fcb2c3](https://github.com/DrudgeCAS/gristmill/commit/5fcb2c3f00dc420a4370acbacd01e4af7943b78f))
* **build:** Use the current scikit-build-core configuration keys ([#29](https://github.com/DrudgeCAS/gristmill/issues/29)) ([c7e2e23](https://github.com/DrudgeCAS/gristmill/commit/c7e2e23a7269add1be42ea68921733e9972b12ee))
* Replace the unsound pruning rather than only reordering it ([dedb226](https://github.com/DrudgeCAS/gristmill/commit/dedb2267f05b3337ebe9c613869a2c2a2637125d))
* Stop a sum with an index-free product crashing the interpreter ([3fd69a9](https://github.com/DrudgeCAS/gristmill/commit/3fd69a9408c7954d7cb11ddab87990d38af1211f))
* Stop a sum with an index-free product crashing the interpreter ([32cbfb2](https://github.com/DrudgeCAS/gristmill/commit/32cbfb201fc7fd40672562005596703184fbf731)), closes [#45](https://github.com/DrudgeCAS/gristmill/issues/45)
* Stop unordered iteration deciding the optimization result ([#39](https://github.com/DrudgeCAS/gristmill/issues/39)) ([0329393](https://github.com/DrudgeCAS/gristmill/commit/0329393d02453aa316249deeae76eee492e5544f)), closes [#31](https://github.com/DrudgeCAS/gristmill/issues/31)
* Stop vertex numbering deciding which constrictions are found ([d6f7b83](https://github.com/DrudgeCAS/gristmill/commit/d6f7b8362c7bb7bd127338a7f56e3cef537bcc84))
* Stop vertex numbering deciding which constrictions are found ([be4f5a1](https://github.com/DrudgeCAS/gristmill/commit/be4f5a19a66023c578f302e727545d30ce7fbb9b)), closes [#43](https://github.com/DrudgeCAS/gristmill/issues/43)
* Take the vertex numbering out of the last two places it was read ([fafd8fe](https://github.com/DrudgeCAS/gristmill/commit/fafd8fe58417f762379c6fe4b2cb4935472725cb))
* Verify a result whatever its external symbols are called ([1d24ddc](https://github.com/DrudgeCAS/gristmill/commit/1d24ddcfe9d905826943dc1f89e16d9171cd6aba)), closes [#49](https://github.com/DrudgeCAS/gristmill/issues/49)


### Performance Improvements

* Bring the search back to master's speed on coupled-cluster equations ([e4f6889](https://github.com/DrudgeCAS/gristmill/commit/e4f688938ff1d7b939c3168af2a5c52af4bdc019))


### Dependencies

* Bump fbitset and libparenth, and cover the crash they fixed ([#36](https://github.com/DrudgeCAS/gristmill/issues/36)) ([70fdb56](https://github.com/DrudgeCAS/gristmill/commit/70fdb566d500f1d6b1880ef26665fe24d3f40d09))
* **deps:** bump deps/fbitset from `765b2f7` to `bf87654` ([1a6432d](https://github.com/DrudgeCAS/gristmill/commit/1a6432d318761379934b5a2805fef0c89f927fee))
* **deps:** bump deps/libparenth from `a1e4899` to `b1e240e` ([a4bb367](https://github.com/DrudgeCAS/gristmill/commit/a4bb367e773f0bbcce5eb5cbb2ad8f813c5064d2))
* **deps:** bump filelock from 3.20.0 to 3.20.3 ([#42](https://github.com/DrudgeCAS/gristmill/issues/42)) ([863b408](https://github.com/DrudgeCAS/gristmill/commit/863b408b2ab3af5b017fbdb85b3207dc42a34f00))
* **deps:** bump idna from 3.10 to 3.15 ([#32](https://github.com/DrudgeCAS/gristmill/issues/32)) ([e35f9ae](https://github.com/DrudgeCAS/gristmill/commit/e35f9ae6b11080a50b0c49b71ac08baca6421fa5))
* **deps:** bump pygments from 2.19.2 to 2.20.0 ([#40](https://github.com/DrudgeCAS/gristmill/issues/40)) ([aa2acc8](https://github.com/DrudgeCAS/gristmill/commit/aa2acc81256dd8faf851db95018ab78c05970ab8))
* **deps:** bump pytest from 9.0.2 to 9.0.3 ([#33](https://github.com/DrudgeCAS/gristmill/issues/33)) ([cd63afb](https://github.com/DrudgeCAS/gristmill/commit/cd63afbb0cc62cfebef79667daaa266782d4485a))
* **deps:** bump requests from 2.32.5 to 2.33.0 ([#41](https://github.com/DrudgeCAS/gristmill/issues/41)) ([4e87f64](https://github.com/DrudgeCAS/gristmill/commit/4e87f64133b591b29b1ee2ff2e04238f68df999e))
* **deps:** bump the dependencies group across 1 directory with 4 updates ([#28](https://github.com/DrudgeCAS/gristmill/issues/28)) ([98d5b13](https://github.com/DrudgeCAS/gristmill/commit/98d5b13ac6da4715c496160ee6ad3e4c95725309))
* **deps:** bump the dependencies group with 6 updates ([#26](https://github.com/DrudgeCAS/gristmill/issues/26)) ([a4f0eb6](https://github.com/DrudgeCAS/gristmill/commit/a4f0eb6510162ccb8ea04ed23daee4053ce9a746))
* **deps:** bump tornado from 6.5.2 to 6.5.7 ([#34](https://github.com/DrudgeCAS/gristmill/issues/34)) ([d02e179](https://github.com/DrudgeCAS/gristmill/commit/d02e1791f4da73fc25df5360423a02ae1b0b4a65))
* **deps:** bump urllib3 from 2.5.0 to 2.7.0 ([#35](https://github.com/DrudgeCAS/gristmill/issues/35)) ([de45ccf](https://github.com/DrudgeCAS/gristmill/commit/de45ccf717e251c8e752f3d2c8e7c019ead25cbb))
* Point cpypp at the merged fix ([7ca0bb1](https://github.com/DrudgeCAS/gristmill/commit/7ca0bb1272e4a12ef7ef28811a033ba8a9627dec))
* Refresh the pinned drudge revision ([#30](https://github.com/DrudgeCAS/gristmill/issues/30)) ([22ab2c7](https://github.com/DrudgeCAS/gristmill/commit/22ab2c7797b3fdf9835ba89491f62ce79f3581cb))


### Code Refactoring

* Build each tie key once, and correct the orientation comment ([c26c240](https://github.com/DrudgeCAS/gristmill/commit/c26c24037033b2cdb6e7c961b4439ea11c43acb2))


### Tests

* Pin that an unrelated term does not disturb the rest ([4713b1d](https://github.com/DrudgeCAS/gristmill/commit/4713b1d3fd74850ce446235060c039e8d95f166e))
* Pin the result against the memoir walk order and keep rand_constr random ([6692241](https://github.com/DrudgeCAS/gristmill/commit/6692241604e3bf834088667356d7989843879797))


### Build System

* Migrate to scikit build system ([#22](https://github.com/DrudgeCAS/gristmill/issues/22)) ([7316e49](https://github.com/DrudgeCAS/gristmill/commit/7316e491debe373400fbd5df0781d23e7d3d81f4))
