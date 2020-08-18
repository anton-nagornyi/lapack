// @ts-ignore
export const offset = <T extends Float32Array | Int32Array>(arr: T, offset: number): T =>
{
    // @ts-ignore
    if (arr instanceof Float32Array) return new Float32Array(arr.buffer, offset * 4);
    // @ts-ignore
    if (arr instanceof Int32Array) return new Int32Array(arr.buffer, offset * 4);
};