import typescript from 'rollup-plugin-typescript2';
import {terser} from "rollup-plugin-terser";
import dts from "rollup-plugin-dts";

export default [
	{
		input: 'src/index.ts',
		output: {
			dir: './',
			format: 'cjs',
		},
		plugins: [typescript({
			include: [
				"./src/**/*.ts"
			],
			useTsconfigDeclarationDir: true,
		}),
			terser()],
	},
	{
		input: "./@types/index.d.ts",
		output: [{ file: "./index.d.ts", format: "es" }],
		plugins: [dts()],
	},
];