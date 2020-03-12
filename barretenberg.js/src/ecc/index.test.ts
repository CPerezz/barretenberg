import { BarretenbergWasm } from '../wasm';

describe('field_muls', () => {
  let barretenberg!: BarretenbergWasm;

  beforeAll(async () => {
    barretenberg = new BarretenbergWasm();
    await barretenberg.init();
  });

  it('time field muls', async () => {
    const start = Date.now();
    barretenberg.exports().do_n_muls(100000000, 0);
    const result = Buffer.from(barretenberg.getMemory().slice(0, 32));
    const duration = Date.now() - start;
    console.log(result);
    console.log('Time: ', duration);
  });
});
