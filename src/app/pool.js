// Worker 池: 任务按"通道"分组 (main 晕面 / obs 焦散 / shell 球壳 / cone 锥切片), 各通道独立取消.
// 同一通道再次 run() 时, 该通道排队的任务丢弃, 正在跑的 Worker 直接终止重建 (即时取消); 其余通道不受影响.
// Worker 源码由构建脚本内联成字符串 (WORKER_SRC), 用 Blob URL 启动 — 单文件 HTML 双击即可离线运行.

export class Pool {
  constructor(src, size) {
    this.url = URL.createObjectURL(new Blob([src], { type: 'text/javascript' }));
    this.size = Math.max(1, size);
    this.workers = [];        // { w, busy: task | null }
    this.queue = [];          // task = { ch, gen, kind, args, resolve, reject }
    this.gen = {};            // 通道 → 当前代号
    for (let i = 0; i < this.size; i++) this._spawn(i);
  }

  _spawn(i) {
    const w = new Worker(this.url);
    const slot = { w, busy: null };
    w.onmessage = e => {
      const t = slot.busy;
      slot.busy = null;
      if (t && t.gen === this.gen[t.ch]) {
        const m = e.data;
        if (m.ok) t.resolve(m.data); else t.reject(new Error(m.error));
      }
      this._pump();
    };
    w.onerror = e => {
      e.preventDefault();
      const t = slot.busy;
      slot.busy = null;
      if (t && t.gen === this.gen[t.ch]) t.reject(new Error(e.message || 'Worker error'));
      this._pump();
    };
    this.workers[i] = slot;
  }

  _pump() {
    for (let i = 0; i < this.workers.length && this.queue.length; i++) {
      const s = this.workers[i];
      if (s.busy) continue;
      const t = this.queue.shift();
      s.busy = t;
      s.w.postMessage({ kind: t.kind, args: t.args });
    }
  }

  /** 取消通道 ch 的全部任务 (排队的丢弃, 在跑的终止). */
  cancel(ch) {
    this.gen[ch] = (this.gen[ch] || 0) + 1;
    this.queue = this.queue.filter(t => t.ch !== ch);
    for (let i = 0; i < this.workers.length; i++) {
      const s = this.workers[i];
      if (s.busy && s.busy.ch === ch) { s.w.terminate(); this._spawn(i); }
    }
  }

  /**
   * 在通道 ch 上跑一批任务 (先取消该通道旧任务). tasks = [{ kind, args }], 排在前面的先发.
   * onEach(i, data) 每完成一个回调一次; 返回 Promise (全部完成 → results 数组; 被取消 → null).
   */
  run(ch, tasks, onEach) {
    this.cancel(ch);
    const gen = this.gen[ch];
    const results = new Array(tasks.length);
    let left = tasks.length;
    return new Promise((resolve, reject) => {
      if (left === 0) { resolve(results); return; }
      tasks.forEach((tk, i) => {
        this.queue.push({
          ch, gen, kind: tk.kind, args: tk.args,
          resolve: d => {
            results[i] = d;
            if (onEach) onEach(i, d);
            if (--left === 0) resolve(results);
          },
          reject: err => { this.cancel(ch); reject(err); },
        });
      });
      this._pump();
    }).then(r => (gen === this.gen[ch] ? r : null));
  }
}
