// 界面语言: 简体中文 / English. 文字在 strings.js ({ 键: [中文, English] }) 与 terms.js (术语解释).
// t('键', {变量}); HTML 里带 data-i18n / data-i18n-title / -aria / -ph 的元素由 applyI18n 换字.

import STR from './strings.js';
import TERMS from './terms.js';

let cur = 'zh';

export function detectLang(saved) {
  if (saved === 'zh' || saved === 'en') return saved;
  const first = String((navigator.languages && navigator.languages[0]) || navigator.language || '').toLowerCase();
  return first.startsWith('zh') ? 'zh' : 'en';
}
export function setLang(code) {
  cur = code === 'en' ? 'en' : 'zh';
  document.documentElement.lang = cur === 'zh' ? 'zh-CN' : 'en';
  return cur;
}
export const getLang = () => cur;

export function t(key, vars) {
  const e = STR[key];
  let s = e === undefined ? key : (cur === 'en' ? e[1] : e[0]);
  if (vars) s = s.replace(/\{(\w+)\}/g, (m, k) => (vars[k] !== undefined ? vars[k] : m));
  return s;
}

/** 术语: { title, body(HTML) } */
export function term(key) {
  const e = TERMS[key];
  if (!e) return null;
  const v = cur === 'en' ? e.en : e.zh;
  return { title: v[0], body: v[1] };
}
export const termKeys = () => Object.keys(TERMS);

export function applyI18n(root = document) {
  for (const el of root.querySelectorAll('[data-i18n]')) {
    const s = t(el.dataset.i18n);
    if (s.includes('<')) el.innerHTML = s; else el.textContent = s;
  }
  for (const [attr, name] of [['data-i18n-title', 'title'], ['data-i18n-aria', 'aria-label'], ['data-i18n-ph', 'placeholder']])
    for (const el of root.querySelectorAll(`[${attr}]`)) el.setAttribute(name, t(el.getAttribute(attr)));
}
