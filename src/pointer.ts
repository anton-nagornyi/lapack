export class Pointer<T>
{
    constructor(value: T)
    {
        this.value = value;
    }
    value: T;

    assign(other: Pointer<T>)
    {
        this.value = other.value;
    }
}